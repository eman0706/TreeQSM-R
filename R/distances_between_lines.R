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

#' Distances Between Lines
#'
#' Distanza minima (perpendicolare comune) tra una retta 3D singola
#' (P1, d1) e una o piu' rette (P2mat, d2mat, righe = rette diverse).
#' Formula standard delle rette sghembe (closest point of approach).
#'
#' @param P1 Punto sulla prima retta (lunghezza 3)
#' @param d1 Direzione della prima retta (lunghezza 3, non deve essere unitaria)
#' @param P2mat Punti sulle altre rette (n x 3)
#' @param d2mat Direzioni delle altre rette (n x 3)
#' @return Lista con Dist (distanza minima), DistOnLine1 (parametro sulla
#'   prima retta), DistOnLines (parametro sulle altre rette - la
#'   posizione del punto piu' vicino, misurata da P2mat lungo d2mat)
#' @export
distances_between_lines <- function(P1, d1, P2mat, d2mat) {
  if (is.null(dim(P2mat))) P2mat <- matrix(P2mat, nrow = 1)
  if (is.null(dim(d2mat))) d2mat <- matrix(d2mat, nrow = 1)
  n <- nrow(P2mat)
  
  d1m <- matrix(d1, n, 3, byrow = TRUE)
  P1m <- matrix(P1, n, 3, byrow = TRUE)
  
  w0 <- P1m - P2mat
  a <- rowSums(d1m * d1m)
  b <- rowSums(d1m * d2mat)
  cc <- rowSums(d2mat * d2mat)
  d <- rowSums(d1m * w0)
  e <- rowSums(d2mat * w0)
  
  denom <- a * cc - b * b
  denom[abs(denom) < 1e-12] <- 1e-12  # rette quasi parallele
  
  t1 <- (b * e - cc * d) / denom
  t2 <- (a * e - b * d) / denom
  
  closest1 <- P1m + t1 * d1m
  closest2 <- P2mat + t2 * d2mat
  Dist <- sqrt(rowSums((closest1 - closest2)^2))
  
  list(Dist = Dist, DistOnLine1 = t1, DistOnLines = t2)
}
