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

#' Matrix-Vector Subtraction
#'
#' Subtracts a row vector from every row of a matrix (row-wise broadcast).
#' R equivalent of TreeQSM's mat_vec_subtraction.m
#' (mat - repmat(vec, size(mat,1), 1)).
#'
#' @param mat n x m matrix
#' @param vec Vector of length m
#' @return n x m matrix, mat with vec subtracted from every row
#' @export
mat_vec_subtraction <- function(mat, vec) {
  mat <- as.matrix(mat)
  mat - matrix(vec, nrow(mat), ncol(mat), byrow = TRUE)
}