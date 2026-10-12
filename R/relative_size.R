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

#' Saturate to uint8
#'
#' Replica il comportamento di MATLAB per una variabile dichiarata
#' 'uint8': ogni assegnazione viene arrotondata e saturata all'intervallo
#' [0,255] (non un clamp esplicito a parte - avviene ad OGNI singola
#' assegnazione nell'originale, dato che RS e' dichiarato uint8 fin
#' dall'inizio).
#'
#' @param x Vettore numerico
#' @return Vettore intero saturato in [0,255]
sat_uint8 <- function(x) {
  x <- round(x)
  x[x < 0]   <- 0
  x[x > 255] <- 255
  as.integer(x)
}

#' Relative Size
#'
#' Porting fedele di relative_size.m (Version 2.00, Copyright (C)
#' 2014-2017 Pasi Raumonen) - sostituisce la versione precedente
#' ("Simplified relative_size", un'invenzione basata solo sull'altezza
#' media del segmento, mai una vera traduzione fedele nonostante fosse
#' classificata come tale in precedenza).
#'
#' Usa la segmentazione esistente e la sua struttura di ramificazione
#' per determinare la dimensione relativa dei cover set distribuiti sul
#' nuovo cover. L'idea e' che la dimensione relativa decresca man mano
#' che la dimensione del ramo decresce: alla base di un ramo e'
#' proporzionale alla dimensione della base del tronco (misurata come
#' numero di cover set nei primi layer), mentre verso la punta del ramo
#' decresce fino al minimo. Dimensione massima teorica 256 alla base
#' del tronco (in pratica saturata a 255, essendo RS un uint8 in
#' MATLAB), minima 1 alla punta di ogni ramo.
#'
#' NOTA DI FEDELTA': l'originale calcola l'ordine di ramificazione con
#' "order = order+order" (RADDOPPIO, non incremento "+1") - con order
#' iniziale 0, questo lo mantiene sempre 0 per costruzione (0+0=0
#' indefinitamente), rendendo la parte "branching order" dell'algoritmo
#' effettivamente un no-op (probabile refuso nell'originale pubblicato,
#' forse voleva essere "order+1"). Fedeli al comportamento REALE del
#' codice pubblicato, non a quello presumibilmente inteso - la sessione
#' ha gia' visto un caso analogo (i pesi in cylinder_fitting.R) dove
#' "correggere" un apparente bug isolato, senza poter ricalibrare tutto
#' il resto del sistema attorno ad esso, ha peggiorato nettamente i
#' risultati.
#'
#' @param P Point cloud completa (Nx3)
#' @param cover Cover set (ball, center, neighbor)
#' @param segment Segmentazione (segments, ChildSegment)
#' @return RS, vettore intero (1-256, saturato a 255) di lunghezza nrow(P)
#' @export
relative_size <- function(P, cover, segment) {
  Bal  <- cover$ball
  Cen  <- cover$center
  Nei  <- cover$neighbor
  Segs <- segment$segments
  SChi <- segment$ChildSegment
  np <- nrow(P)
  ns <- length(Segs)
  
  if (ns == 0) return(rep(128L, np))
  
  ## Usa ordine di ramificazione e altezza come informazione a priori
  Ord <- integer(ns)
  Cvec <- SChi[[1]]
  Cvec <- Cvec[!is.na(Cvec) & Cvec >= 1 & Cvec <= ns]
  order_val <- 0L
  safety <- 0L
  while (length(Cvec) > 0 && safety < ns + 10L) {
    safety <- safety + 1L
    order_val <- order_val + order_val
    Ord[Cvec] <- order_val
    Cvec <- unlist(SChi[Cvec], use.names = FALSE)
    Cvec <- Cvec[!is.na(Cvec) & Cvec >= 1 & Cvec <= ns]
  }
  maxO <- order_val + 1L
  
  ## Determina l'altezza dell'albero
  Cen_coords <- P[Cen, , drop = FALSE]
  Top <- max(Cen_coords[, 3])
  Bot <- min(Cen_coords[, 3])
  H <- Top - Bot
  if (!is.finite(H) || H <= 0) H <- 1  # protezione contro divisione per zero
  
  ## Determina "base size" rispetto alla base del tronco
  BaseSize <- numeric(ns)
  S1 <- Segs[[1]]
  n1 <- length(S1)
  if (n1 >= 2) {
    m1 <- min(6L, n1)
    BaseSize[1] <- mean(sapply(S1[2:m1], length))
  } else {
    BaseSize[1] <- length(S1[[1]])
  }
  if (!is.finite(BaseSize[1]) || BaseSize[1] <= 0) BaseSize[1] <- 1
  
  if (ns >= 2) {
    for (i in 2:ns) {
      Si <- Segs[[i]]
      ni <- length(Si)
      if (ni >= 2) {
        mi <- min(6L, ni)
        BaseSize[i] <- ceiling(mean(sapply(Si[2:mi], length)) / BaseSize[1] * 256)
      } else {
        BaseSize[i] <- length(Si[[1]]) / BaseSize[1] * 256
      }
      bot_i <- min(Cen_coords[Si[[1]], 3])
      h_i <- bot_i - Bot
      BS <- ceiling(256 * (maxO - Ord[i]) / maxO * (H - h_i) / H)
      if (BaseSize[i] > BS) BaseSize[i] <- BS
    }
  }
  BaseSize[1] <- 256
  
  ## Determina la dimensione relativa per i punti (RS e' uint8 in
  ## MATLAB: ogni assegnazione satura implicitamente a [0,255])
  TS <- 1
  RS <- integer(np)
  for (i in seq_len(ns)) {
    Si <- Segs[[i]]
    s <- length(Si)
    if (s == 0) next
    for (j in seq_len(s)) {
      Q <- Si[[j]]
      Q <- Q[!is.na(Q) & Q >= 1 & Q <= length(Bal)]
      if (length(Q) == 0) next
      pts <- unlist(Bal[Q], use.names = FALSE)
      pts <- pts[!is.na(pts) & pts >= 1 & pts <= np]
      if (length(pts) == 0) next
      val <- BaseSize[i] - (BaseSize[i] - TS) * sqrt((j - 1) / s)
      RS[pts] <- sat_uint8(val)
    }
  }
  
  ## Aggiusta la dimensione relativa alla base dei segmenti figli
  RS0 <- RS
  for (i in seq_len(ns)) {
    C <- SChi[[i]]
    C <- C[!is.na(C) & C >= 1 & C <= ns]
    n <- length(C)
    if (n > 0) {
      for (j in seq_len(n)) {
        Sj <- Segs[[C[j]]]
        if (length(Sj) == 0) next
        B <- Sj[[1]]
        B <- B[!is.na(B) & B >= 1 & B <= length(Nei)]
        if (length(B) == 0) next
        N <- unlist(Nei[B], use.names = FALSE)
        if (length(Sj) > 1) {
          N <- setdiff(N, Sj[[2]])
        }
        N <- union(N, B)
        N <- N[!is.na(N) & N >= 1 & N <= length(Bal)]
        if (length(N) == 0) next
        Npts <- unlist(Bal[N], use.names = FALSE)
        Npts <- Npts[!is.na(Npts) & Npts >= 1 & Npts <= np]
        if (length(Npts) == 0) next
        RS[Npts] <- sat_uint8(RS0[Npts] / 2)
      }
    }
  }
  
  RS[RS < 1L] <- 1L
  RS
}
