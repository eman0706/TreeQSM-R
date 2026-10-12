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

#' Segment Direction
#'
#' Porting fedele di segment_direction (sotto-funzione di
#' correct_segments.m). Stima la direzione locale di un segmento
#' attorno al layer "nl", usando layer un po' sopra e un po' sotto
#' (fino a 3 di distanza, saltando quelli vuoti).
#'
#' @param Ce Centri dei cover set
#' @param Seg Un singolo segmento (lista di layer)
#' @param nl Indice del layer attorno a cui stimare la direzione
#' @return Vettore direzione (lunghezza 3), o c(0,0,0) se non determinabile
#' @export
segment_direction <- function(Ce, Seg, nl) {
  bot <- if (nl - 3 > 0) nl - 3 else 1
  j <- 1L
  while (j < 3L && length(Seg[[bot]]) == 0L) {
    bot <- bot + 1L
    j <- j + 1L
  }
  top <- if (nl + 2 <= length(Seg)) nl + 2 else length(Seg)
  j <- 1L
  while (j < 3L && length(Seg[[top]]) == 0L) {
    top <- top - 1L
    j <- j + 1L
  }
  
  if (top > bot) {
    Bot <- average(Ce[Seg[[bot]], , drop = FALSE])
    Top <- average(Ce[Seg[[top]], , drop = FALSE])
    V <- Top - Bot
    D <- V / sqrt(sum(V^2))
  } else {
    D <- c(0, 0, 0)
  }
  D
}