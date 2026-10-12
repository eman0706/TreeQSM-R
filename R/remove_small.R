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

#' Remove Small
#'
#' Porting fedele di remove_small (sotto-funzione di
#' correct_segments.m). Rimuove segmenti figli piccoli/poco chiari: per
#' ciascun figlio, confronta la sua estensione radiale con quella del
#' genitore nella stessa zona - se sono geometricamente indistinguibili
#' (il figlio non si allontana significativamente dall'asse del
#' genitore), lo considera parte del genitore stesso e lo rimuove
#' (ricorsivamente, controllando anche i nipoti).
#'
#' @param Ce Centri dei cover set
#' @param Segs segment$segments
#' @param SPar segment$ParentSegment (matrice ns x 2)
#' @param SChi segment$ChildSegment
#' @return Lista con Segs, SPar, SChi aggiornati
#' @export
remove_small <- function(Ce, Segs, SPar, SChi) {
  Segment <- Segs[[1]]
  ns <- length(Segment)
  EndL <- if (ns > 10) 10L else ns
  End   <- average(Ce[Segment[[EndL]], , drop = FALSE])
  Start <- average(Ce[Segment[[1]], , drop = FALSE])
  V <- End - Start
  V <- V / sqrt(sum(V^2))
  Sets0 <- unlist(Segment[seq_len(EndL)])
  MaxRad <- max(distances_to_line(Ce[Sets0, , drop = FALSE], V, Start)$d)
  
  Nseg <- length(Segs)
  Keep <- rep(TRUE, Nseg)
  
  for (i in seq_len(Nseg)) {
    if (Keep[i]) {
      ChildSegs <- SChi[[i]]
      if (length(ChildSegs) > 0) {
        n <- length(ChildSegs)
        Segment <- Segs[[i]]
        ns <- length(Segment)
        for (j in seq_len(n)) {
          nl <- SPar[ChildSegs[j], 2]
          StartL <- if (nl > 10) nl - 10L else 1L
          EndL   <- if (ns - nl > 10) nl + 10L else ns
          End   <- average(Ce[Segment[[EndL]], , drop = FALSE])
          Start <- average(Ce[Segment[[StartL]], , drop = FALSE])
          V <- End - Start
          V <- V / sqrt(sum(V^2))
          
          ChildSets <- unlist(Segs[[ChildSegs[j]]])
          
          distChild <- max(distances_to_line(Ce[ChildSets, , drop = FALSE], V, Start)$d)
          
          if (distChild < MaxRad + 0.06) {
            ParentSets <- unlist(Segment[StartL:EndL])
            
            distPar <- max(distances_to_line(Ce[ParentSets, , drop = FALSE], V, Start)$d)
            if ((distChild - distPar < 0.02) ||
                (distChild / distPar < 1.2 && distChild - distPar < 0.06)) {
              ChildChildSegs <- SChi[[ChildSegs[j]]]
              nc <- length(ChildChildSegs)
              if (nc == 0) {
                # Rimuovi, nessun segmento figlio-del-figlio
                Keep[ChildSegs[j]] <- FALSE
                Segs[[ChildSegs[j]]] <- list()
                SPar[ChildSegs[j], ] <- c(0L, 0L)
                SChi[[i]] <- setdiff(ChildSegs, ChildSegs[j])
              } else {
                L <- unlist(SChi[ChildChildSegs])
                if (length(L) == 0) {
                  Jflag <- rep(FALSE, nc)
                  for (k in seq_len(nc)) {
                    segment_k <- Segs[[ChildChildSegs[k]]]
                    if (length(segment_k) == 0) {
                      Jflag[k] <- TRUE
                    } else {
                      segment1 <- c(unlist(segment_k), ParentSets)
                      distSeg <- max(distances_to_line(Ce[segment1, , drop = FALSE], V, Start)$d)
                      if ((distSeg - distPar < 0.02) ||
                          (distSeg / distPar < 1.2 && distSeg - distPar < 0.06)) {
                        Jflag[k] <- TRUE
                      }
                    }
                  }
                  if (all(Jflag)) {
                    # Rimuovi il figlio e tutti i suoi figli diretti
                    ChildChildSegs1 <- c(ChildChildSegs, ChildSegs[j])
                    for (cc in ChildChildSegs1) {
                      Segs[[cc]] <- list()
                      SChi[[cc]] <- integer(0)
                    }
                    Keep[ChildChildSegs1] <- FALSE
                    SPar[ChildChildSegs1, ] <- matrix(0L, length(ChildChildSegs1), 2)
                    SChi[[i]] <- setdiff(ChildSegs, ChildSegs[j])
                  }
                }
              }
            }
          }
        }
      }
      if (i == 1L) {
        MaxRad <- MaxRad / 2
      }
    }
  }
  
  # Modifica i segmenti e la loro indicizzazione
  Segs <- Segs[Keep]
  n_keep <- sum(Keep)
  Ind <- seq_len(Nseg)
  Ind[Keep]  <- seq_len(n_keep)
  Ind[!Keep] <- 0L
  SPar <- SPar[Keep, , drop = FALSE]
  Jp <- SPar[, 1] > 0
  if (any(Jp)) SPar[Jp, 1] <- Ind[SPar[Jp, 1]]
  
  SChi_kept <- SChi
  for (i in seq_len(Nseg)) {
    if (Keep[i]) {
      ChildSegs2 <- SChi_kept[[i]]
      if (length(ChildSegs2) > 0) {
        mapped <- Ind[ChildSegs2]
        mapped <- mapped[mapped != 0]
        SChi_kept[[i]] <- mapped
      } else {
        SChi_kept[[i]] <- integer(0)
      }
    }
  }
  SChi <- SChi_kept[Keep]
  
  list(Segs = Segs, SPar = SPar, SChi = SChi)
}