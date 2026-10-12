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

#' Verticalcat
#'
#' Equivalente di verticalcat.m: concatena i cover set di tutti i layer
#' di un segmento in un unico vettore "Sets", registrando in "IndSets"
#' (nl x 2) l'indice di inizio/fine di ciascun layer dentro Sets.
#' @param Seg Lista di layer (ciascuno un vettore di indici di cover set)
#' @return Lista con Sets (vettore concatenato) e IndSets (nl x 2)
#' @export
verticalcat <- function(Seg) {
  nl <- length(Seg)
  lens <- sapply(Seg, length)
  ends <- cumsum(lens)
  starts <- c(1L, if (nl > 1) head(ends, -1) + 1L else integer(0))
  Sets <- unlist(Seg)
  IndSets <- cbind(starts, ends)
  list(Sets = Sets, IndSets = IndSets)
}
