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

#' Search Stem Top
#'
#' Porting fedele di search_stem_top (sotto-funzione di
#' correct_segments.m, Pasi Raumonen). Cerca il segmento che, se fuso
#' con il segmento 1 (tronco) seguendo la catena genitore-figlio,
#' produce il tronco "piu' dritto" possibile: tra i segmenti la cui
#' punta e' vicina al punto piu' alto dell'albero (entro SearchDist,
#' che si allarga se nessuno e' abbastanza vicino), sceglie quello con
#' il minor rapporto lunghezza/distanza-lineare (cioe' il meno curvo).
#'
#' NOTA DI IMPLEMENTAZIONE: il puntatore "a" (posizione di scrittura nel
#' buffer "Segment") e' tracciato esplicitamente con lo stesso identico
#' pattern di incremento dell'originale MATLAB (buffer pre-allocato,
#' scrittura a blocchi, "a" che punta sempre alla PROSSIMA posizione
#' libera) - deliberatamente NON semplificato con un append dinamico
#' ricalcolato al volo, per eliminare ogni rischio di disallineamento
#' nei confronti "SP >= a-2" rispetto all'originale.
#'
#' @param P Point cloud (n x 3)
#' @param Ce Centri dei cover set (matrice, allineata a cover$center)
#' @param Bal cover$ball
#' @param Segs segment$segments
#' @param SPar segment$ParentSegment (matrice ns x 2)
#' @param dmin Diametro minimo di riferimento (per il passo di
#'   campionamento della lunghezza approssimata)
#' @return Indice del segmento che forma la cima del tronco
#' @export
search_stem_top <- function(P, Ce, Bal, Segs, SPar, dmin) {
  nseg <- length(Segs)
  SegHeight <- numeric(nseg)
  HorDist   <- numeric(nseg)
  s <- Segs[[1]][[1]]
  StemCen <- average(Ce[s, , drop = FALSE])
  
  for (i in seq_len(nseg)) {
    last_layer <- Segs[[i]][[length(Segs[[i]])]]
    S <- last_layer[1]
    SegHeight[i] <- Ce[S, 3]
    HorDist[i] <- sqrt(sum((Ce[S, 1:2] - StemCen[1:2])^2))
  }
  
  Top <- max(SegHeight)
  HeiDist <- Top - SegHeight
  Dist <- sqrt(HorDist^2 + HeiDist^2)
  
  LenDisRatio    <- 2
  SearchDist     <- 0.5
  MaxLenDisRatio <- 1.05
  StemTop        <- 1L
  BUFLEN <- 3000L
  safety_outer   <- 0L
  
  while (LenDisRatio > MaxLenDisRatio) {
    StemTops <- seq_len(nseg)
    I <- Dist < SearchDist
    safety_inner <- 0L
    while (!any(I)) {
      SearchDist <- SearchDist + 0.5
      I <- Dist < SearchDist
      safety_inner <- safety_inner + 1L
      if (safety_inner > 10000L) {
        warning("search_stem_top: ciclo di ricerca SearchDist non converge (probabili NaN in Dist) - uso fallback StemTop=1")
        return(1L)
      }
    }
    StemTops <- StemTops[I]
    
    n <- length(StemTops)
    Stems <- vector("list", n)
    
    for (j in seq_len(n)) {
      Seg  <- Segs[[1]]
      spar <- SPar
      
      if (StemTops[j] != 1L) {
        SubSegs <- integer(200)
        SubSegs[1] <- StemTops[j]
        nsegs <- 1L
        segment <- StemTops[j]
        safety_chain <- 0L
        while (segment != 1L) {
          segment <- SPar[segment, 1]
          nsegs <- nsegs + 1L
          SubSegs[nsegs] <- segment
          safety_chain <- safety_chain + 1L
          if (safety_chain > nseg + 10L) {
            warning("search_stem_top: risalita verso la radice non converge (probabile ciclo in SPar) - uso fallback StemTop=1")
            return(1L)
          }
        }
        
        # Buffer pre-allocato, fedele a "Segment = cell(3000,1)"
        Segment <- vector("list", BUFLEN)
        a <- length(Seg)
        Segment[seq_len(a)] <- Seg
        a <- a + 1L
        
        if (nsegs > 2L) {
          for (i2 in seq_len(nsegs - 2L)) {
            I_seg <- SubSegs[nsegs - i2]
            J_seg <- SubSegs[nsegs - i2 - 1L]
            SP <- spar[I_seg, 2]
            SegC <- Segs[[I_seg]]
            sp <- spar[J_seg, 2]
            if (SP >= a - 2L) {
              Segment[a:(a + sp - 1L)] <- SegC[seq_len(sp)]
              spar[J_seg, 2] <- a + sp - 1L
              a <- a + sp
            } else {
              Segment[(SP + 1L):(SP + sp)] <- SegC[seq_len(sp)]
              a <- SP + sp + 1L
              spar[J_seg, 2] <- SP + sp
            }
            SubSegs[nsegs - i2] <- 1L
          }
        }
        
        I_seg <- SubSegs[1]
        SP <- spar[I_seg, 2]
        SegC <- Segs[[I_seg]]
        nc <- length(SegC)
        if (SP >= a - 2L) {
          Segment[a:(a + nc - 1L)] <- SegC
          a <- a + nc - 1L
        } else {
          Segment[(SP + 1L):(SP + nc)] <- SegC
          a <- SP + nc
        }
        Stems[[j]] <- Segment[seq_len(a)]
      } else {
        Stems[[j]] <- Seg
      }
    }
    
    N <- ceiling(0.5 / dmin / 1.4)
    Lengths <- numeric(n)
    Heights <- numeric(n)
    
    for (i2 in seq_len(n)) {
      Seg2 <- Stems[[i2]]
      ns2  <- length(Seg2)
      m <- if (ceiling(ns2 / N) > floor(ns2 / N)) ceiling(ns2 / N) else ceiling(ns2 / N) + 1L
      Nodes <- matrix(0, m, 3)
      for (j2 in seq_len(m)) {
        I_idx <- (j2 - 1L) * N + 1L
        if (I_idx > ns2) I_idx <- ns2
        S2 <- Seg2[[I_idx]]
        if (length(S2) > 1L) {
          Nodes[j2, ] <- average(Ce[S2, , drop = FALSE])
        } else {
          S2b <- Bal[[S2]]
          Nodes[j2, ] <- average(P[S2b, , drop = FALSE])
        }
      }
      if (m > 1) {
        V <- Nodes[2:m, , drop = FALSE] - Nodes[1:(m - 1), , drop = FALSE]
        Lengths[i2] <- sum(sqrt(rowSums(V * V)))
        V2 <- Nodes[m, ] - Nodes[1, ]
        Heights[i2] <- sqrt(sum(V2 * V2))
      } else {
        Lengths[i2] <- 0
        Heights[i2] <- 0
      }
    }
    
    ratio <- Lengths / Heights
    if (length(ratio) == 0) {
      warning("search_stem_top: nessun candidato StemTops valido in questo giro - uso fallback StemTop corrente")
      break
    }
    I_min <- which.min(ratio)
    if (length(I_min) == 0) {
      warning("search_stem_top: which.min(ratio) vuoto (probabile ratio tutto NaN) - uso fallback StemTop corrente")
      break
    }
    LenDisRatio <- ratio[I_min]
    if (!is.finite(LenDisRatio)) {
      warning("search_stem_top: LenDisRatio non finito (probabile Heights=0) - uso fallback StemTop corrente")
      break
    }
    StemTop <- StemTops[I_min]
    SearchDist <- SearchDist + 1
    if (SearchDist > 3) {
      MaxLenDisRatio <- 1.1
      if (SearchDist > 5) {
        MaxLenDisRatio <- 1.15
        if (SearchDist > 7) {
          MaxLenDisRatio <- 5
        }
      }
    }
    safety_outer <- safety_outer + 1L
    if (safety_outer > 100L) {
      warning("search_stem_top: ciclo esterno non converge dopo 100 iterazioni - uso il candidato corrente")
      break
    }
  }
  
  StemTop
}