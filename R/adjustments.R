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

#' Adjustments
#'
#' Porting fedele COMPLETO di adjustments.m (parte di cylinders.m
#' v3.0.0, Pasi Raumonen), verificato riga per riga contro il vero
#' sorgente MATLAB (recuperato integralmente in una sessione
#' successiva a quella in cui questo file era stato dichiarato
#' "parziale"). Applica correzioni finali ai cilindri di un segmento:
#' tetto/pavimento del raggio basato sul cilindro genitore, correzione
#' di rastremazione (taper - lineare per rami molto sottili o corti,
#' parabola ai minimi quadrati per il caso generale), ottimizzazione
#' dei punti di partenza (fit circolare 2D), continuita' tra cilindri
#' consecutivi, e collegamento del primo cilindro al genitore se
#' troppo lontano. Nessuna divergenza trovata rispetto all'originale -
#' solo protezioni difensive aggiuntive (tryCatch attorno a solve(),
#' controllo di norma positiva prima di dividere).
#'
#' @param cyl Cilindri del segmento corrente (radius, length, start,
#'   axis, radius0, SurfCov, mad - vettori/matrici allineati)
#' @param parcyl Cilindro genitore (radius, length, start, axis - liste
#'   di lunghezza 1, o liste vuote se non c'e' genitore)
#' @param inputs Lista di input (usa MinCylRad, ParentCor, TaperCor)
#' @param Regs Lista delle regioni di punti (matrici Nx3) usate per
#'   ciascun cilindro fittato - output $Reg di cylinder_fitting()
#' @return cyl aggiornato
#' @export
adjustments <- function(cyl, parcyl, inputs, Regs) {
  
  nc <- length(cyl$radius)
  if (nc == 0) return(cyl)
  Mod <- rep(FALSE, nc)
  SC  <- cyl$SurfCov
  
  MinCylRad <- inputs$MinCylRad %||% 0.005
  
  ## ---- Determina il raggio massimo e minimo ----
  has_parent <- !is.null(parcyl$radius) && length(parcyl$radius) > 0 &&
    is.finite(parcyl$radius[1])
  if (has_parent) {
    MaxR <- 0.95 * parcyl$radius[1]
    MaxR <- max(MaxR, MinCylRad)
  } else {
    a <- min(3, nc)
    MaxR <- 1.25 * max(cyl$radius[seq_len(a)])
  }
  
  MinR <- if (any(SC > 0.7)) min(cyl$radius[SC > 0.7]) else NA_real_
  if (!is.na(MinR) && min(cyl$radius) < MinR / 2) {
    MinR <- if (any(SC > 0.4)) min(cyl$radius[SC > 0.4]) else NA_real_
  } else if (is.na(MinR)) {
    MinR <- if (any(SC > 0.4)) min(cyl$radius[SC > 0.4]) else MinCylRad
  }
  if (is.na(MinR)) MinR <- MinCylRad
  
  ## ---- Applica i tetti massimo/minimo ----
  I <- cyl$radius < MinR
  cyl$radius[I] <- MinR
  Mod[I] <- TRUE
  
  ParentCor <- isTRUE(inputs$ParentCor %||% TRUE)
  if (ParentCor || nc <= 3) {
    I <- (cyl$radius > MaxR & SC < 0.7) | (cyl$radius > 1.2 * MaxR)
    cyl$radius[I] <- MaxR
    Mod[I] <- TRUE
    if (nc <= 3) {
      I <- (cyl$radius > 0.75 * MaxR & SC < 0.7)
      if (any(I)) {
        r <- pmax(SC[I] / 0.7 * cyl$radius[I], MinR)
        cyl$radius[I] <- r
        Mod[I] <- TRUE
      }
    }
  }
  
  ## ---- Correzione locale: bassa copertura tra due cilindri ben coperti ----
  if (nc > 2) {
    for (i in 2:(nc - 1)) {
      if (SC[i] < 0.7 && SC[i - 1] >= 0.7 && SC[i + 1] >= 0.7) {
        cyl$radius[i] <- 0.5 * (cyl$radius[i - 1] + cyl$radius[i + 1])
        Mod[i] <- TRUE
      }
    }
  }
  
  ## ---- Correzione di rastremazione (taper) ----
  TaperCor <- isTRUE(inputs$TaperCor %||% TRUE)
  if (TaperCor) {
    if (max(cyl$radius) < 0.001) {
      ## Rami molto sottili: rastremazione lineare forzata
      if (nc > 2) {
        r <- sort(cyl$radius)
        r <- r[2:(length(r) - 1)]
        a <- 2 * mean(r)
        if (a > max(r)) a <- min(0.01, max(r))
        b <- min(0.5 * min(cyl$radius), 0.001)
        cyl$radius <- seq(a, b, length.out = nc)
      } else if (nc > 1) {
        r <- max(cyl$radius)
        cyl$radius <- c(r, 0.5 * r)
      }
      Mod <- rep(TRUE, nc)
      
    } else if (nc > 4) {
      ## Parabola di rastremazione (fit ai minimi quadrati)
      branchlen <- sum(cyl$length[seq_len(nc)])
      L <- cyl$length / 2 + c(0, cumsum(cyl$length[seq_len(nc - 1)]))
      Taper1 <- c(L, branchlen)
      Taper2 <- c(1.05 * cyl$radius, MinR)
      sc <- c(SC, 1)
      
      A_mat <- matrix(c(
        sum(sc * Taper1^4), sum(sc * Taper1^2),
        sum(sc * Taper1^2), sum(sc)
      ), 2, 2, byrow = TRUE)
      y_vec <- c(sum(sc * Taper2 * Taper1^2), sum(sc * Taper2))
      x_sol <- tryCatch(solve(A_mat, y_vec), error = function(e) c(-0.0001, MinR))
      x_sol[1] <- min(x_sol[1], -0.0001)
      
      Ru <- x_sol[1] * L^2 + x_sol[2]
      Ru[Ru < MinR] <- MinR
      if (max(Ru) > MaxR) {
        a <- max(Ru)
        Ru <- MaxR / a * Ru
      }
      Rl <- 0.75 * Ru
      Rl[Rl < MinR] <- MinR
      
      I <- cyl$radius > Ru & SC < 0.7
      cyl$radius[I] <- Ru[I] + (cyl$radius[I] - Ru[I]) * SC[I] / 0.7
      Mod[I] <- TRUE
      I <- cyl$radius > 1.333 * Ru & SC >= 0.7
      cyl$radius[I] <- Ru[I] + (cyl$radius[I] - Ru[I]) * SC[I]
      Mod[I] <- TRUE
      I <- (cyl$radius < Rl & SC < 0.7) | (cyl$radius < 0.5 * Rl)
      cyl$radius[I] <- Rl[I]
      Mod[I] <- TRUE
      
    } else {
      ## Rami corti: rastremazione lineare tra estremi ben coperti
      R <- cyl$radius
      if (sum(SC >= 0.7) > 1) {
        a <- max(R[SC >= 0.7]); b <- min(R[SC >= 0.7])
      } else if (sum(SC >= 0.7) == 1) {
        a <- max(R[SC >= 0.7]); b <- min(R)
      } else {
        a <- sum(R * SC / sum(SC)); b <- min(R)
      }
      Ru <- seq(a, b, length.out = nc)
      I <- SC < 0.7 & !Mod
      cyl$radius[I] <- Ru[I] + (R[I] - Ru[I]) * SC[I] / 0.7
      Mod[I] <- TRUE
    }
  }
  
  ## ---- Modifica i punti di partenza ottimizzandoli per il raggio e
  ## l'asse dati (least_squares_circle_centre + surface_coverage2) ----
  nr <- length(Regs)
  Reg <- NULL
  for (i in seq_len(nc)) {
    if (Mod[i]) {
      if (nr == nc) {
        Reg <- Regs[[i]]
      } else if (i > 1) {
        Reg <- Regs[[i - 1]]
      }
      if (!is.null(Reg) && abs(cyl$radius[i] - cyl$radius0[i]) > 0.005 &&
          (nr == nc || (nr < nc && i > 1))) {
        P <- sweep(Reg, 2, cyl$start[i, ])
        ow <- orthonormal_vectors(cyl$axis[i, ])
        U <- ow$V; V <- ow$W
        P_proj <- cbind(as.vector(P %*% U), as.vector(P %*% V))
        cir <- tryCatch(
          least_squares_circle_centre(P_proj, c(0, 0), cyl$radius[i]),
          error = function(e) NULL
        )
        if (!is.null(cir) && isTRUE(cir$conv) && isTRUE(cir$rel)) {
          cyl$start[i, ] <- cyl$start[i, ] + cir$point[1] * U + cir$point[2] * V
          cyl$mad[i] <- cir$mad
          dl <- distances_to_line(Reg, cyl$axis[i, ], cyl$start[i, ])
          h  <- dl$h
          if (min(h) < -0.001) {
            cyl$length[i] <- max(h) - min(h)
            cyl$start[i, ] <- cyl$start[i, ] + min(h) * cyl$axis[i, ]
            dl <- distances_to_line(Reg, cyl$axis[i, ], cyl$start[i, ])
            h  <- dl$h
          }
          a  <- max(0.02, 0.2 * cyl$radius[i])
          nl <- ceiling(cyl$length[i] / a)
          nl <- max(nl, 4)
          ns <- ceiling(2 * pi * cyl$radius[i] / a)
          ns <- max(ns, 10)
          ns <- min(ns, 36)
          cyl$SurfCov[i] <- surface_coverage2(cyl$axis[i, ], cyl$length[i], dl$V, h, nl, ns)
        }
      }
    }
  }
  
  ## ---- Continuita': ogni cilindro deve iniziare dove finisce il
  ## precedente ----
  if (nc > 1) {
    for (j in 2:nc) {
      U <- cyl$start[j, ] - cyl$start[j - 1, ] - cyl$length[j - 1] * cyl$axis[j - 1, ]
      if (sqrt(sum(U^2)) > 0.0001) {
        N <- cyl$axis[j, ]
        if (sqrt(sum(N^2)) > 0) {
          ow <- orthonormal_vectors(N)
          Vv <- ow$V; Ww <- ow$W
          M <- cbind(N, Vv, Ww)
          x <- tryCatch(solve(M, U), error = function(e) c(0, 0, 0))
          cyl$start[j, ] <- cyl$start[j, ] - x[1] * N
          if (x[1] > 0) {
            cyl$length[j] <- cyl$length[j] + x[1]
          } else if (cyl$length[j] + x[1] > 0) {
            cyl$length[j] <- cyl$length[j] + x[1]
          }
        }
      }
    }
  }
  
  ## ---- Collega il primo cilindro al genitore se troppo lontano ----
  if (has_parent) {
    dl <- distances_to_line(matrix(cyl$start[1, ], 1, 3), parcyl$axis[1, ], parcyl$start[1, ])
    d <- dl$d[1] - parcyl$radius[1]
    if (is.finite(d) && d > 0.001) {
      Vn <- dl$V[1, ]
      Vn <- parcyl$radius[1] * Vn / sqrt(sum(Vn^2))
      h <- dl$h[1]
      B <- dl$B[1, ]
      taper <- cyl$start[1, ]
      E <- taper + cyl$length[1] * cyl$axis[1, ]
      if (h >= 0 && h <= parcyl$length[1]) {
        cyl$start[1, ] <- parcyl$start[1, ] + B + Vn
      } else if (h < 0) {
        cyl$start[1, ] <- parcyl$start[1, ] + Vn
      } else {
        cyl$start[1, ] <- parcyl$start[1, ] + parcyl$length[1] * parcyl$axis[1, ] + Vn
      }
      new_axis <- E - cyl$start[1, ]
      new_len  <- sqrt(sum(new_axis^2))
      if (new_len > 1e-9) {
        cyl$axis[1, ]  <- new_axis / new_len
        cyl$length[1]  <- new_len
      }
    }
  }
  
  cyl
}