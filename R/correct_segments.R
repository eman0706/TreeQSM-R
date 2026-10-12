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

#' Correct Segments
#'
#' Porting fedele della vera correct_segments.m (Pasi Raumonen).
#' Sostituisce la vecchia patch "fatta in casa" mai verificata contro il
#' sorgente originale - probabile causa originale della frammentazione
#' del tronco su piu' segmenti inseguita per tutta la sessione.
#'
#' Rende il tronco e i rami "il piu' lunghi possibile" (modify_topology,
#' la parte piu' importante - fonde segmenti che sono la continuazione
#' fisica dello stesso ramo), poi opzionalmente rimuove segmenti piccoli
#' e poco chiari (remove_small), poi opzionalmente espande le basi dei
#' rami nel genitore (modify_parent). Infine costruisce
#' SegmentOfPoint/branch1indexes/branch2indexes/branch3indexes - i campi
#' che define_main_branches() usera' per riparare la connettivita' del
#' cover nel secondo passaggio.
#'
#' NOTA DI FEDELTA': l'originale trasforma ParentSegment da matrice
#' (ns x 2) a vettore semplice (solo la prima colonna) in output - qui
#' preservato cosi' com'e', dato che serve in questo formato a
#' define_main_branches(). Il secondo passaggio di segments() ricostruisce
#' comunque la propria ParentSegment (matrice) da zero, quindi non c'e'
#' conflitto con il resto della pipeline.
#'
#' @param P Point cloud
#' @param cover Cover set (ball, center, neighbor)
#' @param segment Segmenti (da segments())
#' @param inputs Parametri (usa PatchDiam1, PatchDiam2Max)
#' @param RemSmall Se TRUE, rimuove segmenti piccoli poco chiari
#' @param ModBases Se TRUE, modifica le basi dei segmenti
#' @param AddChild Se TRUE (e ModBases=TRUE), la base espansa viene
#'   aggiunta al figlio; se FALSE, viene solo rimossa dal genitore
#' @return segment aggiornato, con i nuovi campi SegmentOfPoint,
#'   branch1indexes, branch2indexes, branch3indexes
#' @export
correct_segments <- function(P, cover, segment, inputs,
                             RemSmall = TRUE, ModBases = FALSE, AddChild = FALSE) {
  
  Bal  <- cover$ball
  Segs <- segment$segments
  SPar <- segment$ParentSegment
  SChi <- segment$ChildSegment
  Ce   <- P[cover$center, , drop = FALSE]
  
  if (is.null(dim(SPar))) SPar <- cbind(as.integer(SPar), 0L)
  
  ## Rende tronco e rami il piu' lunghi possibile
  if (isTRUE(RemSmall)) {
    mt <- modify_topology(P, Ce, Bal, Segs, SPar, SChi, inputs$PatchDiam2Max)
  } else {
    mt <- modify_topology(P, Ce, Bal, Segs, SPar, SChi, inputs$PatchDiam1)
  }
  Segs <- mt$Segs; SPar <- mt$SPar; SChi <- mt$SChi
  
  ## Rimuove i segmenti figli piccoli e poco chiari
  if (isTRUE(RemSmall)) {
    rs <- remove_small(Ce, Segs, SPar, SChi)
    Segs <- rs$Segs; SPar <- rs$SPar; SChi <- rs$SChi
  }
  
  # Controlla la coerenza dei vettori vuoti
  ns <- length(Segs)
  for (i in seq_len(ns)) {
    if (length(SChi[[i]]) == 0) SChi[[i]] <- integer(0)
  }
  
  if (isTRUE(ModBases)) {
    ns <- length(Segs)
    # FIX (trappola MATLAB "for i=2:ns"): con ns==1, "2:1" in R darebbe
    # c(2,1) invece del range vuoto che MATLAB produrrebbe - guardia
    # esplicita per sicurezza.
    if (ns >= 2) {
      if (isTRUE(AddChild)) {
        # Aggiunge la base espansa al figlio, la rimuove dal genitore
        for (i in 2:ns) {
          SegC <- Segs[[i]]
          SegP <- Segs[[SPar[i, 1]]]
          mp <- modify_parent(P, Bal, Ce, SegP, SegC, SPar[i, 2], inputs$PatchDiam1)
          Segs[[SPar[i, 1]]] <- mp$SegP
          SegC[[1]] <- mp$Base
          Segs[[i]] <- SegC
        }
      } else {
        # Rimuove solo la base espansa dal genitore
        for (i in 2:ns) {
          if (SPar[i, 1] > 1) {
            SegC <- Segs[[i]]
            SegP <- Segs[[SPar[i, 1]]]
            mp <- modify_parent(P, Bal, Ce, SegP, SegC, SPar[i, 2], inputs$PatchDiam2Max)
            Segs[[SPar[i, 1]]] <- mp$SegP
          }
        }
      }
    }
  }
  
  # NOTA FEDELTA': ParentSegment diventa vettore semplice qui (solo la
  # prima colonna), come nell'originale.
  SPar_vec <- SPar[, 1]
  
  segment$segments      <- Segs
  segment$ParentSegment <- SPar_vec
  segment$ChildSegment  <- SChi
  
  ## Genera i dati di segmento per i punti
  np <- nrow(P)
  ns <- length(Segs)
  SegmentOfPoint <- integer(np)
  for (i in seq_len(ns)) {
    S <- Segs[[i]]
    all_cs <- unlist(S)
    all_cs <- all_cs[!is.na(all_cs) & all_cs >= 1L & all_cs <= length(Bal)]
    if (length(all_cs) == 0) next
    pts <- unlist(Bal[all_cs])
    pts <- pts[!is.na(pts) & pts >= 1L & pts <= np]
    SegmentOfPoint[pts] <- i
  }
  segment$SegmentOfPoint <- SegmentOfPoint
  
  C <- SChi[[1]]
  segment$branch1indexes <- C
  if (length(C) > 0) {
    C <- unlist(SChi[C])
    segment$branch2indexes <- C
    if (length(C) > 0) {
      C <- unlist(SChi[C])
      segment$branch3indexes <- C
    } else {
      segment$branch3indexes <- integer(0)
    }
  } else {
    segment$branch2indexes <- integer(0)
    segment$branch3indexes <- integer(0)
  }
  
  segment
}