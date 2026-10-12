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

#' Distances to a Line
#'
#' Computes, for each point in Q, the perpendicular (shortest) distance to
#' an infinite 3D line defined by a point on the line and a direction
#' vector, together with the signed position of the projected point along
#' the line direction (measured from LinePoint).
#'
#' @param Q Points to measure (n x 3), or a single point (length 3)
#' @param LineDirec Direction vector of the line (length 3, need not be
#'   pre-normalized)
#' @param LinePoint A point on the line (length 3)
#' @return List with:
#'   \item{d}{Perpendicular distance from each point in Q to the line}
#'   \item{h}{Signed distance, along LineDirec, from LinePoint to the
#'             projection of each point onto the line}
#'   \item{V}{The perpendicular vector(s) from the line to each point (n x 3)}
#' @export
distances_to_line <- function(Q, LineDirec, LinePoint) {
  single_point <- is.null(dim(Q))
  Qm <- if (single_point) matrix(Q, nrow = 1) else as.matrix(Q)
  n <- nrow(Qm)
  
  dir_len <- norm_vec(LineDirec)
  if (!is.finite(dir_len) || dir_len < 1e-12) {
    stop("distances_to_line: LineDirec has (near) zero length")
  }
  u <- as.numeric(LineDirec) / dir_len
  
  W <- Qm - matrix(LinePoint, n, 3, byrow = TRUE)
  h <- as.numeric(W %*% u)
  Proj <- matrix(h, n, 3) * matrix(u, n, 3, byrow = TRUE)
  V <- W - Proj
  d <- norm_vec(V)
  
  list(d = d, h = h, V = V)
}
