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

#' Surface Coverage 2
#'
#' Porting fedele di surface_coverage2.m (Pasi Raumonen). Calcola la
#' copertura di superficie (numero tra 0 e 1) dei punti sulla
#' superficie di un cilindro definito da Axis/Len - versione piu'
#' semplice/veloce di surface_coverage_filtering.R, SENZA filtro dei
#' punti (qui si vuole solo il numero, non un sottoinsieme filtrato).
#'
#' @param Axis Asse del cilindro (vettore lunghezza 3)
#' @param Len Lunghezza del cilindro
#' @param Vec Vettori dai punti all'asse (matrice Nx3)
#' @param height Altezze dei punti dalla base del cilindro (vettore N)
#' @param nl Numero di layer (divisioni lungo l'asse)
#' @param ns Numero di settori (divisioni angolari)
#' @return SurfCov, numero tra 0 e 1
#' @export
surface_coverage2 <- function(Axis, Len, Vec, height, nl, ns) {
  ow <- orthonormal_vectors(Axis)
  U <- ow$V
  W <- ow$W
  
  VecUW <- cbind(as.vector(Vec %*% U), as.vector(Vec %*% W))
  ang <- atan2(VecUW[, 2], VecUW[, 1]) + pi
  
  I <- ceiling(height / Len * nl)
  I[I == 0] <- 1
  I[I > nl] <- nl
  
  Jc <- ceiling(ang / (2 * pi) * ns)
  Jc[Jc == 0] <- 1
  
  K <- I + (Jc - 1) * nl
  
  length(unique(K)) / nl / ns
}
