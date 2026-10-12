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

#' Point Model Distance
#'
#' Porting fedele di point_model_distance.m (Version 2.1.1, Copyright (C)
#' 2015-2021 Pasi Raumonen). Calcola le distanze dei punti reali alla
#' superficie del modello a cilindri - usata per la metrica di
#' ottimizzazione (default 'all_mean_dis' in select_optimum) e per
#' valutare la qualita' complessiva del fit.
#'
#' NOTA DI FEDELTA': aggiunti controlli difensivi (clamp degli indici
#' del cubo di ricerca entro i limiti della griglia, skip se il range
#' risulta vuoto) non espliciti nell'originale - l'originale si affida
#' al margine NE per garantire indici sempre validi; qui aggiungiamo
#' una protezione esplicita per evitare errori di indicizzazione in R
#' in casi limite, senza cambiare il comportamento nei casi normali.
#'
#' @param P Point cloud completa (Nx3)
#' @param cylinder Output di cylinders() (radius, length, start, axis,
#'   BranchOrder)
#' @return Lista pmdistance con CylDist, mean/median/max/std (globale,
#'   Trunk, Branch, Branch1, Branch2)
#' @export
point_model_distance <- function(P, cylinder) {
  Rad  <- cylinder$radius
  Len  <- cylinder$length
  Sta  <- cylinder$start
  Axe  <- cylinder$axis
  BOrd <- cylinder$BranchOrder
  
  ## Seleziona casualmente il 25% (o max 1 milione) dei punti
  np0 <- nrow(P)
  a   <- min(0.25 * np0, 1000000)
  I_samp <- as.logical(round(0.5 / (1 - a / np0) * runif(np0)))
  P <- P[I_samp, , drop = FALSE]
  
  ## Partiziona i punti in cubi
  L  <- 2 * median(Len)
  NE <- max(3, min(10, ceiling(max(Len) / L))) + 3
  cp <- cubical_partition(P, L, NE)
  Partition <- cp$Partition
  Info <- cp$Info
  Min  <- Info[1:3]
  EL   <- Info[7]
  NE   <- Info[8]
  
  ## Coordinate cubiche dei punti di partenza dei cilindri
  CC <- cbind(
    floor((Sta[, 1] - Min[1]) / EL),
    floor((Sta[, 2] - Min[2]) / EL),
    floor((Sta[, 3] - Min[3]) / EL)
  ) + NE + 1
  
  ## Numero di cubi necessari per ciascun cilindro
  Ncubes <- ceiling(Len / L)
  
  ## Corregge N affinche' gli indici dei cubi non siano troppo piccoli o grandi
  Ibad <- CC[, 1] < Ncubes + 1
  Ncubes[Ibad] <- CC[Ibad, 1] - 1
  Ibad <- CC[, 2] < Ncubes + 1
  Ncubes[Ibad] <- CC[Ibad, 2] - 1
  Ibad <- CC[, 3] < Ncubes + 1
  Ncubes[Ibad] <- CC[Ibad, 3] - 1
  Ibad <- CC[, 1] + Ncubes + 1 > Info[4]
  Ncubes[Ibad] <- Info[4] - CC[Ibad, 1] - 1
  Ibad <- CC[, 2] + Ncubes + 1 > Info[5]
  Ncubes[Ibad] <- Info[5] - CC[Ibad, 2] - 1
  Ibad <- CC[, 3] + Ncubes + 1 > Info[6]
  Ncubes[Ibad] <- Info[6] - CC[Ibad, 3] - 1
  Ncubes <- pmax(Ncubes, 0)
  
  ## Calcola le distanze ai cilindri
  n  <- length(Rad)
  np <- nrow(P)
  Dist <- matrix(0, np, 2)  # colonna 1: distanza, colonna 2: cilindro piu' vicino
  Dist[, 1] <- 2  # distanza grande iniziale
  Data <- vector("list", n)
  N_dim <- dim(Partition)
  
  for (i in seq_len(n)) {
    x1 <- max(1, CC[i, 1] - Ncubes[i]); x2 <- min(CC[i, 1] + Ncubes[i], N_dim[1])
    y1 <- max(1, CC[i, 2] - Ncubes[i]); y2 <- min(CC[i, 2] + Ncubes[i], N_dim[2])
    z1 <- max(1, CC[i, 3] - Ncubes[i]); z2 <- min(CC[i, 3] + Ncubes[i], N_dim[3])
    if (x1 > x2 || y1 > y2 || z1 > z2) next
    Par <- Partition[x1:x2, y1:y2, z1:z2]
    points <- unlist(Par, use.names = FALSE)
    if (length(points) == 0) next
    
    dl <- distances_to_line(P[points, , drop = FALSE], Axe[i, ], Sta[i, ])
    d  <- abs(dl$d - Rad[i])
    h  <- dl$h
    Data[[i]] <- cbind(d, h, points)
    
    I_close <- d < Dist[points, 1]
    J <- h >= 0
    K <- h <= Len[i]
    Lg <- d < 0.5
    M <- I_close & J & K & Lg
    pts_upd <- points[M]
    if (length(pts_upd) > 0) {
      Dist[pts_upd, 1] <- d[M]
      Dist[pts_upd, 2] <- i
    }
  }
  
  ## Calcola le distanze per i punti non ancora assegnati (non "sul
  ## lato" del cilindro, ma vicino alle estremita')
  for (i in seq_len(n)) {
    if (!is.null(Data[[i]])) {
      d <- Data[[i]][, 1]
      h <- Data[[i]][, 2]
      points <- Data[[i]][, 3]
      I_close <- d < Dist[points, 1]
      J <- h >= -0.1 & h <= 0
      K <- h <= Len[i] + 0.1 & h >= Len[i]
      Lg <- d < 0.5
      M <- I_close & (J | K) & Lg
      pts_upd <- points[M]
      if (length(pts_upd) > 0) {
        Dist[pts_upd, 1] <- d[M]
        Dist[pts_upd, 2] <- i
      }
    }
  }
  
  ## Seleziona solo il 95% piu' vicino delle distanze per ciascun cilindro
  assigned <- Dist[, 2]
  Cyl <- split(seq_len(np), factor(assigned, levels = seq_len(n)))
  DistCyl <- numeric(n)
  for (i in seq_len(n)) {
    idx <- Cyl[[i]]
    m <- length(idx)
    if (m > 19) {
      d_sorted <- sort(Dist[idx, 1])
      DistCyl[i] <- mean(d_sorted[seq_len(floor(0.95 * m))])
    } else if (m > 0) {
      DistCyl[i] <- mean(Dist[idx, 1])
    }
  }
  
  ## Costruisce l'output
  pmdistance <- list()
  pmdistance$CylDist <- DistCyl
  pmdistance$median  <- median(DistCyl)
  pmdistance$mean    <- mean(DistCyl)
  pmdistance$max     <- max(DistCyl)
  pmdistance$std     <- sd(DistCyl)
  
  Ttrunk <- BOrd == 0
  B1 <- BOrd == 1
  B2 <- BOrd == 2
  Bvals  <- DistCyl[!Ttrunk]
  Tvals  <- DistCyl[Ttrunk]
  B1vals <- DistCyl[B1]
  B2vals <- DistCyl[B2]
  
  pmdistance$TrunkMedian <- median(Tvals)
  pmdistance$TrunkMean   <- mean(Tvals)
  pmdistance$TrunkMax    <- max(Tvals)
  pmdistance$TrunkStd    <- sd(Tvals)
  
  if (length(Bvals) > 0) {
    pmdistance$BranchMedian <- median(Bvals)
    pmdistance$BranchMean   <- mean(Bvals)
    pmdistance$BranchMax    <- max(Bvals)
    pmdistance$BranchStd    <- sd(Bvals)
  } else {
    pmdistance$BranchMedian <- 0
    pmdistance$BranchMean   <- 0
    pmdistance$BranchMax    <- 0
    pmdistance$BranchStd    <- 0
  }
  
  if (length(B1vals) > 0) {
    pmdistance$Branch1Median <- median(B1vals)
    pmdistance$Branch1Mean   <- mean(B1vals)
    pmdistance$Branch1Max    <- max(B1vals)
    pmdistance$Branch1Std    <- sd(B1vals)
  } else {
    pmdistance$Branch1Median <- 0
    pmdistance$Branch1Mean   <- 0
    pmdistance$Branch1Max    <- 0
    pmdistance$Branch1Std    <- 0
  }
  
  if (length(B2vals) > 0) {
    pmdistance$Branch2Median <- median(B2vals)
    pmdistance$Branch2Mean   <- mean(B2vals)
    pmdistance$Branch2Max    <- max(B2vals)
    pmdistance$Branch2Std    <- sd(B2vals)
  } else {
    pmdistance$Branch2Median <- 0
    pmdistance$Branch2Mean   <- 0
    pmdistance$Branch2Max    <- 0
    pmdistance$Branch2Std    <- 0
  }
  
  pmdistance
}
