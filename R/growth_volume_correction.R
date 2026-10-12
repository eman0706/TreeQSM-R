# This file is part of TREEQSM.
# 
# TREEQSM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# TREEQSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

#' Allometry Function
#'
#' Modello allometrico usato da growth_volume_correction: Rad = a*GV^b+c
#' @param x Vettore di parametri c(a, b, c)
#' @param xdata Vettore di growth volume
#' @return Raggio predetto
allometry_fn <- function(x, xdata) {
  x[1] * xdata^x[2] + x[3]
}

#' Growth Volume Correction
#'
#' Porting fedele di growth_volume_correction.m (Version 2.0.0,
#' Copyright (C) 2013-2021 Pasi Raumonen). Usa l'approccio allometrico
#' "growth volume" (= il volume totale "supportato" da un cilindro,
#' introdotto da Jan Hackenberg nei paper "SimpleTree") per correggere
#' il raggio di cilindri troppo grandi o troppo piccoli. Correzione
#' OPZIONALE, disattivata di default (inputs$GrowthVolCor=FALSE) - non
#' fa parte del flusso principale a meno che l'utente non la attivi
#' esplicitamente.
#'
#' NOTA DI FEDELTA': l'originale usa lsqcurvefit (trust-region-
#' reflective, MATLAB Optimization Toolbox) per il fit allometrico
#' Rad=a*GV^b+c. Qui si usa nls() di R (Gauss-Newton) con gli stessi
#' valori iniziali [0.5, 0.5, 0] - non e' bit-esattamente identico
#' all'algoritmo MATLAB, ma converge a valori molto simili per un
#' problema di fitting standard come questo.
#'
#' @param cylinder Output di cylinders() (radius, length, parent,
#'   extension)
#' @param inputs Lista di input (usa GrowthVolFac)
#' @return cylinder con radius eventualmente corretto
#' @export
growth_volume_correction <- function(cylinder, inputs) {
  Rad  <- as.numeric(cylinder$radius)
  Rad0 <- Rad
  Len  <- as.numeric(cylinder$length)
  CPar <- cylinder$parent
  CExt <- cylinder$extension
  n    <- length(Rad)
  
  if (n == 0) return(cylinder)
  
  cat("----------\n")
  cat("Growth volume based correction of cylinder radii:\n")
  
  initial_volume <- round(1000 * pi * sum(Rad^2 * Len))
  cat(sprintf(" Initial_volume (L): %d\n", initial_volume))
  
  ## Definisce i cilindri figli per ciascun cilindro
  ind  <- seq_len(n)
  CChi <- vector("list", n)
  for (i in seq_len(n)) {
    CChi[[i]] <- ind[CPar == i]
  }
  
  ## Calcola il growth volume (propagazione bottom-up, dalle foglie
  ## verso la radice)
  GrowthVol <- numeric(n)
  S <- sapply(CChi, length)
  modify <- S == 0
  GrowthVol[modify] <- pi * Rad[modify]^2 * Len[modify]
  parents <- sort(unique(CPar[modify]))
  if (length(parents) > 0 && parents[1] == 0) {
    parents <- parents[-1]
  }
  
  safety <- 0L
  while (length(parents) > 0 && safety < n + 10L) {
    safety <- safety + 1L
    V <- pi * Rad[parents]^2 * Len[parents]
    m <- length(parents)
    for (i in seq_len(m)) {
      GrowthVol[parents[i]] <- V[i] + sum(GrowthVol[CChi[[parents[i]]]])
    }
    parents <- sort(unique(CPar[parents]))
    if (length(parents) > 0 && parents[1] == 0) {
      parents <- parents[-1]
    }
  }
  if (safety >= n + 10L) {
    warning("growth_volume_correction: propagazione growth volume non converge - interrotta con lo stato corrente")
  }
  
  ## Fit dell'allometria: Rad = a*GV^b + c
  fit_data <- data.frame(GV = GrowthVol, R = Rad)
  fit <- tryCatch(
    nls(R ~ a * GV^b + c, data = fit_data,
        start = list(a = 0.5, b = 0.5, c = 0),
        control = nls.control(maxiter = 200, warnOnly = TRUE)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    cat(" Fit allometria fallito - correzione growth volume saltata\n")
    cat("----------\n")
    return(cylinder)
  }
  X <- coef(fit)
  cat(" Allometry model parameters R = a*GV^b+c:\n")
  cat(sprintf("   Multiplier a: %s\n", format(X[1])))
  cat(sprintf("   Exponent b: %s\n", format(X[2])))
  cat(sprintf("   Intersect c: %s\n", format(X[3])))
  
  ## Raggio predetto dall'allometria
  PredRad <- allometry_fn(X, GrowthVol)
  
  ## Corregge i raggi in base alle predizioni: se il raggio supera
  ## fac-volte o e' sotto 1/fac-volte il raggio predetto, correggilo.
  ## Non aumentare mai il raggio alle punte dei rami (CExt==0).
  fac <- inputs$GrowthVolFac
  modify2 <- Rad < PredRad / fac | Rad > fac * PredRad
  modify2[Rad < PredRad / fac & CExt == 0] <- FALSE
  CorRad <- PredRad[modify2]
  
  D <- 0
  if (any(modify2)) {
    R_mod <- Rad[modify2]
    D <- max(abs(R_mod - CorRad))
  }
  
  Rad[modify2] <- CorRad
  cylinder$radius <- Rad
  
  cat(sprintf(" Modified %d of the %d cylinders\n", sum(modify2), n))
  cat(sprintf(" Largest radius change (cm): %s\n", round(1000 * D) / 10))
  corrected_volume <- round(1000 * pi * sum(Rad^2 * Len))
  cat(sprintf(" Corrected volume (L): %d\n", corrected_volume))
  cat(sprintf(" Change in volume (L): %d\n", corrected_volume - initial_volume))
  cat("----------\n")
  
  cylinder
}