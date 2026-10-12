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

#' Partitions point cloud into cubical cells
#'
#' @param P Point cloud, (n_points x 3)-matrix
#' @param EL Length of the cube edges
#' @param NE Number of empty edge layers (default 3)
#' @return List with Partition, CubeCoord, Info
cubical_partition <- function(P, EL, NE = 3) {
  np <- nrow(P)
  
  Min <- apply(P, 2, min)
  Max <- apply(P, 2, max)
  N <- ceiling((Max - Min) / EL) + 2 * NE + 1
  
  if (any(N > 10000)) {
    stop("Too many cubes requested.")
  }
  
  Info <- c(Min, N, EL, NE)
  
  CubeCoord <- matrix(0, np, 3)
  for (i in 1:3) {
    CubeCoord[, i] <- floor((P[, i] - Min[i]) / EL) + NE + 1
    CubeCoord[, i] <- pmax(pmin(CubeCoord[, i], N[i]), 1)
  }
  
  LexOrd <- (CubeCoord[, 1] - 1) + (CubeCoord[, 2] - 1) * N[1] + 
            (CubeCoord[, 3] - 1) * N[1] * N[2]
  SortOrd <- order(LexOrd)
  LexOrd <- LexOrd[SortOrd]
  
  Partition <- array(list(), dim = N)
  
  p <- 1
  while (p <= np) {
    t <- 1
    while ((p + t <= np) && (LexOrd[p] == LexOrd[p + t])) {
      t <- t + 1
    }
    q <- SortOrd[p]
    cx <- CubeCoord[q, 1]
    cy <- CubeCoord[q, 2]
    cz <- CubeCoord[q, 3]
    Partition[cx, cy, cz] <- list(SortOrd[p:(p + t - 1)])
    p <- p + t
  }
  
  return(list(
    Partition = Partition,
    CubeCoord = CubeCoord,
    Info = Info
  ))
}

#' Partitions point cloud into cubical cells (with Cubes output)
#'
#' @param P Point cloud, (n_points x 3)-matrix  
#' @param EL Length of the cube edges
#' @param NE Number of empty edge layers (default 3)
#' @return List with Partition, CubeCoord, Info, and Cubes
cubical_partition_with_cubes <- function(P, EL, NE = 3) {
  np <- nrow(P)
  
  Min <- apply(P, 2, min)
  Max <- apply(P, 2, max)
  N <- ceiling((Max - Min) / EL) + 2 * NE + 1
  
  if (any(N > 10000)) {
    stop("Too many cubes requested.")
  }
  
  Info <- c(Min, N, EL, NE)
  
  CubeCoord <- matrix(0, np, 3)
  for (i in 1:3) {
    CubeCoord[, i] <- floor((P[, i] - Min[i]) / EL) + NE + 1
    CubeCoord[, i] <- pmax(pmin(CubeCoord[, i], N[i]), 1)
  }
  
  LexOrd <- (CubeCoord[, 1] - 1) + (CubeCoord[, 2] - 1) * N[1] + 
            (CubeCoord[, 3] - 1) * N[1] * N[2]
  SortOrd <- order(LexOrd)
  LexOrd <- LexOrd[SortOrd]
  
  Partition <- array(list(), dim = N)
  Cubes <- array(0, dim = N)
  
  p <- 1
  while (p <= np) {
    t <- 1
    while ((p + t <= np) && (LexOrd[p] == LexOrd[p + t])) {
      t <- t + 1
    }
    q <- SortOrd[p]
    cx <- CubeCoord[q, 1]
    cy <- CubeCoord[q, 2]
    cz <- CubeCoord[q, 3]
    Partition[cx, cy, cz] <- list(SortOrd[p:(p + t - 1)])
    Cubes[cx, cy, cz] <- Cubes[cx, cy, cz] + 1
    p <- p + t
  }
  
  return(list(
    Partition = Partition,
    CubeCoord = CubeCoord,
    Info = Info,
    Cubes = Cubes
  ))
}
