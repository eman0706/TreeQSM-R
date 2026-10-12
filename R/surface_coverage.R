# This file is part of the R port of TREEQSM.
# Original TREEQSM: Raumonen et al. (2013), distributed under the GNU General Public License, version 3 or (at your option) any later version.
# R port: Barca, Greco, Tartarino (2026); developed with an AI assistant (see README).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Computes point surface coverage measure
#'
#' @param P Point cloud (n x 3 matrix)
#' @param Axis Axis direction (1 x 3)
#' @param Point Starting point of the cylinder (1 x 3)
#' @param nl Number of layers in the axis direction
#' @param ns Number of sectors
#' @param Dmin Optional minimum point distance from axis
#' @param Dmax Optional maximum point distance from axis
#' @param full Se FALSE (default, come nargout = 1 in MATLAB) calcola e restituisce solo lo scalare SurfCov;
#'   se TRUE restituisce list(SurfCov, Dis, dis, CylVol) (piu' lento).
#' @return Scalare SurfCov, oppure (full = TRUE) lista con SurfCov, Dis, CylVol, dis
surface_coverage <- function(P, Axis, Point, nl, ns, Dmin = NULL, Dmax = NULL, full = FALSE) {
  # Compute the distances and heights of the points
  result <- distances_to_line(P, Axis, Point)
  d <- result$d
  V <- result$V
  h <- result$h
  h <- h - min(h)
  Len <- max(h)
  
  # Filter out points based on the distance to the axis
  if (!is.null(Dmin) || !is.null(Dmax)) {
    Keep <- rep(TRUE, length(d))
    if (!is.null(Dmin)) {
      Keep <- Keep & (d > Dmin)
    }
    if (!is.null(Dmax)) {
      Keep <- Keep & (d < Dmax)
    }
    V <- V[Keep, , drop = FALSE]
    h <- h[Keep]
    d <- d[Keep]
  }
  
  # Compute SurfCov from 4 different baseline directions
  V0 <- V
  # First planar axes
  orthonormal_result <- orthonormal_vectors(Axis)
  U <- orthonormal_result$V   # FIX: orthonormal_vectors() restituisce $V/$W, non $U/$W (con $U si otteneva NULL e V[, 2] andava fuori limite)
  W <- orthonormal_result$W
  # Rotation matrix to rotate the axes
  R <- rotation_matrix(Axis, 2 * pi / ns / 4)
  
  SurfCov <- numeric(4)
  for (i in 1:4) {
    # Rotate the axes
    if (i > 1) {
      U <- R %*% U
      W <- R %*% W
    }
    
    # Compute the angles (sectors) of the points
    V <- V0 %*% cbind(U, W)
    ang <- atan2(V[, 2], V[, 1]) + pi
    
    # Compute lexicographic order (sector,layer) of every point
    Layer <- ceiling(h / Len * nl)
    Layer[Layer <= 0] <- 1
    Layer[Layer > nl] <- nl
    Sector <- ceiling(ang / 2 / pi * ns)
    Sector[Sector <= 0] <- 1
    LexOrd <- Layer + (Sector - 1) * nl
    
    # Compute SurfCov
    Cov <- matrix(0, nl, ns)
    Cov[LexOrd] <- 1
    SurfCov[i] <- sum(Cov != 0) / (nl * ns)
  }
  SurfCov <- max(SurfCov)
  
  # Compute volume estimate
  output <- list(SurfCov = SurfCov)
  
  if (full) {  # nargout > 1
    # Sort according to increasing lexicographic order
    LexOrd_sorted <- sort(LexOrd, index.return = TRUE)
    LexOrd <- LexOrd_sorted$x
    SortOrd <- LexOrd_sorted$ix
    d <- d[SortOrd]
    
    # Compute mean distance of the sector-layer intersections
    Dis <- matrix(0, nl, ns)  # mean distances
    np <- length(LexOrd)  # number of points
    p <- 1
    while (p <= np) {
      t <- 1
      while ((p + t <= np) && (LexOrd[p] == LexOrd[p + t])) {
        t <- t + 1
      }
      Dis[LexOrd[p]] <- mean(d[p:(p + t - 1)])
      p <- p + t
    }
    
    output$Dis <- Dis
    
    if (TRUE) {  # nargout > 2
      # Interpolate missing distances
      D <- Dis
      dis <- Dis
      Dinv <- D[rev(1:nl), ]
      D <- rbind(cbind(Dinv, Dinv, Dinv), cbind(D, D, D), cbind(Dinv, Dinv, Dinv))   # FIX: tiling 3 x 3 come in MATLAB
      Zero <- Dis == 0
      RadMean <- mean(Dis[Dis > 0])
      for (i in 1:nl) {
        for (j in 1:ns) {
          if (Zero[i, j]) {
            d_local <- D[(i + nl - 1):(i + nl + 1), (j + ns - 1):(j + ns + 1)]
            if (sum(d_local > 0) > 1) {
              dis[i, j] <- mean(d_local[d_local > 0])
            } else {
              d_local <- D[(i + nl - 2):(i + nl + 2), (j + ns - 2):(j + ns + 2)]
              if (sum(d_local > 0) > 1) {
                dis[i, j] <- mean(d_local[d_local > 0])
              } else {
                d_local <- D[(i + nl - 3):(i + nl + 3), (j + ns - 3):(j + ns + 3)]
                if (sum(d_local > 0) > 1) {
                  dis[i, j] <- mean(d_local[d_local > 0])
                } else {
                  dis[i, j] <- RadMean
                }
              }
            }
          }
        }
      }
      # Compute the volume estimate
      r <- as.vector(dis)
      CylVol <- 1000 * pi * sum(r^2) / ns * Len / nl
      output$dis <- dis
      output$CylVol <- CylVol
    }
  }
  
  if (full) return(output)
  return(SurfCov)
}