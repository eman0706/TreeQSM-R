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

#' Cylinder Distribution
#'
#' Porting fedele di cylinder_distribution (sotto-funzione di
#' tree_data.m, Pasi Raumonen). Volume, area e lunghezza delle parti
#' legnose in funzione del diametro, zenith o azimuth dei cilindri
#' (classi di 1cm per diametro, 10 gradi per angoli).
#'
#' @param treedata Lista treedata da aggiornare
#' @param cyl Oggetto cylinder (radius, length, axis)
#' @param dist Uno tra "Dia", "Zen", "Azi"
#' @return treedata aggiornato con i campi VolCyl<dist>, AreCyl<dist>, LenCyl<dist>
#' @export
cylinder_distribution <- function(treedata, cyl, dist) {
  if (dist == "Dia") {
    Par <- cyl$radius
    n <- ceiling(max(200 * cyl$radius, na.rm = TRUE))
    a <- 0.005
  } else if (dist == "Zen") {
    Par <- 180 / pi * acos(pmax(-1, pmin(1, cyl$axis[, 3])))
    n <- 18
    a <- 10
  } else if (dist == "Azi") {
    Par <- 180 / pi * atan2(cyl$axis[, 2], cyl$axis[, 1]) + 180
    n <- 36
    a <- 10
  } else {
    stop("cylinder_distribution: dist non riconosciuto")
  }
  if (!is.finite(n) || n < 1) n <- 1
  
  CylDist <- matrix(0, 3, n)
  for (i in seq_len(n)) {
    K <- !is.na(Par) & Par >= (i - 1) * a & Par < i * a
    CylDist[1, i] <- 1000 * pi * sum(cyl$radius[K]^2 * cyl$length[K])
    CylDist[2, i] <- 2 * pi * sum(cyl$radius[K] * cyl$length[K])
    CylDist[3, i] <- sum(cyl$length[K])
  }
  treedata[[paste0("VolCyl", dist)]] <- CylDist[1, ]
  treedata[[paste0("AreCyl", dist)]] <- CylDist[2, ]
  treedata[[paste0("LenCyl", dist)]] <- CylDist[3, ]
  treedata
}

#' Cylinder Height Distribution
#'
#' Porting fedele di cylinder_height_distribution (sotto-funzione di
#' tree_data.m). Volume, area e lunghezza dei cilindri in funzione
#' dell'altezza (classi di 1 metro), con interpolazione lineare per i
#' cilindri che attraversano il confine tra due classi.
#'
#' @param treedata Lista treedata (deve gia' contenere TreeHeight)
#' @param cylinder Oggetto cylinder (radius, length, start, axis)
#' @param ind Vettore di indici 1:nc
#' @return treedata aggiornato con VolCylHei, AreCylHei, LenCylHei
#' @export
cylinder_height_distribution <- function(treedata, cylinder, ind) {
  Rad <- cylinder$radius
  Len <- cylinder$length
  Axe <- cylinder$axis
  
  MaxHei <- ceiling(treedata$TreeHeight)
  if (!is.finite(MaxHei) || MaxHei < 1) MaxHei <- 1
  treedata$VolCylHei <- numeric(MaxHei)
  treedata$AreCylHei <- numeric(MaxHei)
  treedata$LenCylHei <- numeric(MaxHei)
  
  End <- cylinder$start + cbind(Len * Axe[, 1], Len * Axe[, 2], Len * Axe[, 3])
  bot <- min(cylinder$start[, 3])
  B <- cylinder$start[, 3] - bot
  Top <- End[, 3] - bot
  
  for (j in seq_len(MaxHei)) {
    I1 <- B >= (j - 2) & B < (j - 1)
    J1 <- B >= (j - 1) & B < j
    K1 <- B >= j & B < (j + 1)
    I2 <- Top >= (j - 2) & Top < (j - 1)
    J2 <- Top >= (j - 1) & Top < j
    K2 <- Top >= j & Top < (j + 1)
    C1 <- ind[J1 & J2]
    C2 <- ind[J1 & K2]
    C3 <- ind[J1 & I2]
    C4 <- ind[I1 & J2]
    C5 <- ind[K1 & J2]
    
    v1 <- 1000 * pi * sum(Rad[C1]^2 * Len[C1]); a1 <- 2 * pi * sum(Rad[C1] * Len[C1]); l1 <- sum(Len[C1])
    
    r2 <- (j - B[C2]) / (Top[C2] - B[C2])
    v2 <- 1000 * pi * sum(Rad[C2]^2 * Len[C2] * r2); a2 <- 2 * pi * sum(Rad[C2] * Len[C2] * r2); l2 <- sum(Len[C2] * r2)
    
    r3 <- (B[C3] - j + 1) / (B[C3] - Top[C3])
    v3 <- 1000 * pi * sum(Rad[C3]^2 * Len[C3] * r3); a3 <- 2 * pi * sum(Rad[C3] * Len[C3] * r3); l3 <- sum(Len[C3] * r3)
    
    r4 <- (Top[C4] - j + 1) / (Top[C4] - B[C4])
    v4 <- 1000 * pi * sum(Rad[C4]^2 * Len[C4] * r4); a4 <- 2 * pi * sum(Rad[C4] * Len[C4] * r4); l4 <- sum(Len[C4] * r4)
    
    r5 <- (j - Top[C5]) / (B[C5] - Top[C5])
    v5 <- 1000 * pi * sum(Rad[C5]^2 * Len[C5] * r5); a5 <- 2 * pi * sum(Rad[C5] * Len[C5] * r5); l5 <- sum(Len[C5] * r5)
    
    treedata$VolCylHei[j] <- v1 + v2 + v3 + v4 + v5
    treedata$AreCylHei[j] <- a1 + a2 + a3 + a4 + a5
    treedata$LenCylHei[j] <- l1 + l2 + l3 + l4 + l5
  }
  treedata
}

#' Branch Distribution
#'
#' Porting fedele di branch_distribution (sotto-funzione di
#' tree_data.m). Volume, area, lunghezza e numero di rami in funzione
#' del diametro, altezza, angolo, zenith o azimuth del ramo. Esclude
#' sempre il primo elemento (branch.order[1]), che corrisponde al
#' tronco, non a un vero ramo.
#'
#' @param treedata Lista treedata (deve gia' contenere TreeHeight per dist="Hei")
#' @param branch Oggetto branch (output di branches())
#' @param dist Uno tra "Dia", "Hei", "Ang", "Zen", "Azi"
#' @return treedata aggiornato con VolBranch<dist>, VolBranch1<dist>,
#'   AreBranch<dist>, AreBranch1<dist>, LenBranch<dist>, LenBranch1<dist>,
#'   NumBranch<dist>, NumBranch1<dist>
#' @export
branch_distribution <- function(treedata, branch, dist) {
  n_b <- length(branch$order)
  empty_out <- function(td) {
    for (nm in c("VolBranch", "VolBranch1", "AreBranch", "AreBranch1",
                 "LenBranch", "LenBranch1", "NumBranch", "NumBranch1")) {
      td[[paste0(nm, dist)]] <- numeric(0)
    }
    td
  }
  if (n_b < 2) return(empty_out(treedata))
  
  BOrd <- branch$order[2:n_b]
  BVol <- branch$volume[2:n_b]
  BAre <- branch$area[2:n_b]
  BLen <- branch$length[2:n_b]
  
  if (dist == "Dia") {
    Par <- branch$diameter[2:n_b]
    n <- ceiling(max(100 * Par, na.rm = TRUE))
    a <- 0.005
  } else if (dist == "Hei") {
    Par <- branch$height[2:n_b]
    n <- ceiling(treedata$TreeHeight)
    a <- 1
  } else if (dist == "Ang") {
    Par <- branch$angle[2:n_b]
    n <- 18
    a <- 10
  } else if (dist == "Zen") {
    Par <- branch$zenith[2:n_b]
    n <- 18
    a <- 10
  } else if (dist == "Azi") {
    Par <- branch$azimuth[2:n_b] + 180
    n <- 36
    a <- 10
  } else {
    stop("branch_distribution: dist non riconosciuto")
  }
  if (!is.finite(n) || is.na(n)) n <- 0
  n <- max(0, n)
  
  BranchDist <- matrix(0, 8, max(n, 1))
  if (n > 0) {
    for (i in seq_len(n)) {
      I <- !is.na(Par) & Par >= (i - 1) * a & Par < i * a
      I1 <- I & !is.na(BOrd) & BOrd == 1
      BranchDist[1, i] <- sum(BVol[I])
      BranchDist[2, i] <- sum(BVol[I1])
      BranchDist[3, i] <- sum(BAre[I])
      BranchDist[4, i] <- sum(BAre[I1])
      BranchDist[5, i] <- sum(BLen[I])
      BranchDist[6, i] <- sum(BLen[I1])
      BranchDist[7, i] <- sum(I)
      BranchDist[8, i] <- sum(I1)
    }
    BranchDist <- BranchDist[, seq_len(n), drop = FALSE]
  } else {
    BranchDist <- matrix(0, 8, 0)
  }
  treedata[[paste0("VolBranch", dist)]]  <- BranchDist[1, ]
  treedata[[paste0("VolBranch1", dist)]] <- BranchDist[2, ]
  treedata[[paste0("AreBranch", dist)]]  <- BranchDist[3, ]
  treedata[[paste0("AreBranch1", dist)]] <- BranchDist[4, ]
  treedata[[paste0("LenBranch", dist)]]  <- BranchDist[5, ]
  treedata[[paste0("LenBranch1", dist)]] <- BranchDist[6, ]
  treedata[[paste0("NumBranch", dist)]]  <- BranchDist[7, ]
  treedata[[paste0("NumBranch1", dist)]] <- BranchDist[8, ]
  treedata
}

#' Branch Order Distribution
#'
#' Porting fedele di branch_order_distribution (sotto-funzione di
#' tree_data.m). Volume, area, lunghezza e numero di rami in funzione
#' dell'ordine di ramificazione.
#'
#' @param treedata Lista treedata da aggiornare
#' @param branch Oggetto branch (output di branches())
#' @return treedata aggiornato con VolBranchOrd, AreBranchOrd,
#'   LenBranchOrd, NumBranchOrd
#' @export
branch_order_distribution <- function(treedata, branch) {
  BO <- suppressWarnings(max(branch$order, na.rm = TRUE))
  if (!is.finite(BO)) BO <- 0
  n_rows <- max(1, BO)
  BranchOrdDist <- matrix(0, n_rows, 4)
  for (i in seq_len(n_rows)) {
    I <- !is.na(branch$order) & branch$order == i
    BranchOrdDist[i, 1] <- sum(branch$volume[I])
    BranchOrdDist[i, 2] <- sum(branch$area[I])
    BranchOrdDist[i, 3] <- sum(branch$length[I])
    BranchOrdDist[i, 4] <- sum(I)
  }
  treedata$VolBranchOrd <- BranchOrdDist[, 1]
  treedata$AreBranchOrd <- BranchOrdDist[, 2]
  treedata$LenBranchOrd <- BranchOrdDist[, 3]
  treedata$NumBranchOrd <- BranchOrdDist[, 4]
  treedata
}
