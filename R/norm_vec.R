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

#' Euclidean Norm
#'
#' Computes the Euclidean (L2) norm of a vector. Vectorized: accepts either
#' a single length-3 (or length-n) vector (returns a scalar) or an n x m
#' matrix (returns the row-wise norm as a length-n vector).
#'
#' @param v Vector, or n x m matrix (norm computed per row)
#' @return Norm(s) of v
#' @export
norm_vec <- function(v) {
  if (is.null(dim(v))) {
    sqrt(sum(v^2))
  } else {
    sqrt(rowSums(v^2))
  }
}
