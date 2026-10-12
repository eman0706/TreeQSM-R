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

#' Cross Product
#'
#' Computes the 3D cross product A x B. Vectorized: accepts either two
#' length-3 vectors (returns a length-3 vector) or two n x 3 matrices
#' (returns the row-wise cross product as an n x 3 matrix).
#'
#' @param A First vector(s) (3, or n x 3)
#' @param B Second vector(s) (3, or n x 3)
#' @return Cross product(s) of A and B, same shape convention as input
#' @export
cross_product <- function(A, B) {
  if (is.null(dim(A)) && is.null(dim(B))) {
    # Simple vector case
    c(
      A[2] * B[3] - A[3] * B[2],
      A[3] * B[1] - A[1] * B[3],
      A[1] * B[2] - A[2] * B[1]
    )
  } else {
    # Row-wise matrix case
    A <- matrix(A, ncol = 3)
    B <- matrix(B, ncol = 3)
    if (nrow(A) == 1 && nrow(B) > 1) A <- matrix(A, nrow(B), 3, byrow = TRUE)
    if (nrow(B) == 1 && nrow(A) > 1) B <- matrix(B, nrow(A), 3, byrow = TRUE)
    cbind(
      A[, 2] * B[, 3] - A[, 3] * B[, 2],
      A[, 3] * B[, 1] - A[, 1] * B[, 3],
      A[, 1] * B[, 2] - A[, 2] * B[, 1]
    )
  }
}
