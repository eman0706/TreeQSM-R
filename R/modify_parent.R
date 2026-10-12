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

#' Modify Parent
#'
#' Porting fedele di modify_parent (sotto-funzione di
#' correct_segments.m). Espande la base di un ramo figlio all'indietro
#' nel segmento genitore (per rimuovere "gradini" nel genitore dovuti
#' all'attacco del ramo), poi rimuove l'espansione dal genitore.
#'
#' NOTA DI IMPLEMENTAZIONE: l'originale riceve "base" (un buffer
#' scratch pre-allocato dal chiamante, riusato tra le chiamate) come
#' argomento. Qui e' creato localmente ad ogni chiamata invece che
#' passato dall'esterno - comportamento osservabile identico, dato che
#' l'originale legge sempre e solo la porzione fresca scritta in
#' questa stessa chiamata (base{1:Nlayer+1}), mai dati residui di
#' chiamate precedenti.
#'
#' @param P Point cloud
#' @param Bal cover$ball
#' @param Ce Centri dei cover set
#' @param SegP Segmento genitore (lista di layer)
#' @param SegC Segmento figlio (lista di layer)
#' @param nl Indice del layer nel genitore dove il figlio si attacca
#' @param PatchDiam Diametro di riferimento (PatchDiam1 o PatchDiam2Max
#'   a seconda del chiamante)
#' @return Lista con SegP (genitore aggiornato) e Base (nuova base,
#'   eventualmente espansa, del figlio)
#' @export
modify_parent <- function(P, Bal, Ce, SegP, SegC, nl, PatchDiam) {
  Base <- SegC[[1]]
  
  if (length(Base) > 0) {
    DirChi <- segment_direction(Ce, SegC, 1)
    DirPar <- segment_direction(Ce, SegP, nl)
    
    if (length(Base) > 1) {
      BaseCent <- average(Ce[Base, , drop = FALSE])
      db <- distances_to_line(Ce[Base, , drop = FALSE], DirChi, BaseCent)$d
      DiamBase <- 2 * max(db)
    } else if (length(Bal[[Base]]) > 1) {
      BaseCent <- average(P[Bal[[Base]], , drop = FALSE])
      db <- distances_to_line(P[Bal[[Base]], , drop = FALSE], DirChi, BaseCent)$d
      DiamBase <- 2 * max(db)
    } else {
      BaseCent <- Ce[Base, ]
      DiamBase <- 0
    }
    
    Angle <- abs(sum(DirChi * DirPar))
    Nlayer <- max(3, ceiling(Angle * 2 * DiamBase / PatchDiam))
    if (Nlayer > nl) Nlayer <- nl
    
    layer <- 0L
    base <- vector("list", Nlayer + 2L)
    base[[1]] <- Base
    
    while (layer < Nlayer) {
      Sets <- SegP[[nl - layer]]
      Seg_center <- average(Ce[Sets, , drop = FALSE])
      
      VBase <- Ce[Sets, , drop = FALSE] - matrix(BaseCent, length(Sets), 3, byrow = TRUE)
      h <- as.vector(VBase %*% DirChi)
      B <- matrix(DirChi, length(Sets), 3, byrow = TRUE) * h
      V <- VBase - B
      distSets <- sqrt(rowSums(V * V))
      
      VSeg <- Ce[Sets, , drop = FALSE] - matrix(Seg_center, length(Sets), 3, byrow = TRUE)
      lenBase <- sqrt(rowSums(VBase * VBase))
      lenSeg  <- sqrt(rowSums(VSeg * VSeg))
      
      if (Angle < 0.9) {
        K  <- lenBase < 1.1 / (1 - 0.5 * Angle^2) * lenSeg
        Jm <- distSets < 1.25 * DiamBase
        I  <- K & Jm
      } else {
        I <- distSets < 1.25 * DiamBase
      }
      
      if (all(I) || !any(I)) {
        layer <- Nlayer
      } else {
        SegP[[nl - layer]]   <- Sets[!I]
        base[[layer + 2L]]   <- Sets[I]
        layer <- layer + 1L
      }
    }
    Base <- unlist(base[seq_len(Nlayer + 1L)])
  }
  
  list(SegP = SegP, Base = Base)
}