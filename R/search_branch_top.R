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

#' Search Branch Top
#'
#' Porting fedele di search_branch_top (sotto-funzione di
#' correct_segments.m, Pasi Raumonen). Cerca, tra tutti i sotto-segmenti
#' discendenti dal segmento BI, quello che forma la punta del ramo
#' "piu' dritto": ordina i candidati per distanza lineare (dalla base
#' di BI alla loro punta) decrescente, e sceglie il primo la cui
#' lunghezza approssimata e' vicina alla distanza lineare (ramo poco
#' curvo). Se nessuno soddisfa il criterio stretto, riprova con un
#' criterio piu' permissivo; in ultima istanza usa il candidato con la
#' distanza lineare maggiore.
#'
#' NOTA DI IMPLEMENTAZIONE: stesso approccio di search_stem_top.R - il
#' puntatore "a" (posizione di scrittura nel buffer "Segment") e'
#' tracciato esplicitamente con lo stesso identico pattern di incremento
#' dell'originale MATLAB, non ricalcolato al volo.
#'
#' @param P Point cloud (n x 3)
#' @param Ce Centri dei cover set
#' @param Bal cover$ball
#' @param Segs segment$segments
#' @param SPar segment$ParentSegment (matrice ns x 2)
#' @param SChi segment$ChildSegment
#' @param dmin Diametro minimo di riferimento
#' @param BI Indice del segmento (ramo) di cui cercare la vera punta
#' @return Indice del segmento che forma la punta del ramo
#' @export
search_branch_top <- function(P, Ce, Bal, Segs, SPar, SChi, dmin, BI) {
  ns <- length(Segs)
  BUFLEN <- 3000L
  
  # Definisce tutti i sotto-segmenti di BI (BI stesso + intera discendenza)
  Segments_buf <- integer(ns)
  Segments_buf[1] <- BI
  t <- 2L
  C <- SChi[[BI]]
  n_last <- 0L
  while (length(C) > 0) {
    n_last <- length(C)
    Segments_buf[t:(t + n_last - 1L)] <- C
    C <- unlist(SChi[C])
    t <- t + n_last
  }
  if (t > 2L) {
    t <- t - n_last
  }
  Segments <- Segments_buf[seq_len(t)]
  
  # Distanze lineari dalle punte dei candidati alla base di BI
  LinearDist <- numeric(t)
  Seg0 <- Segs[[Segments[1]]]
  BranchBase <- average(Ce[Seg0[[1]], , drop = FALSE])
  for (i in seq_len(t)) {
    Seg_i <- Segs[[Segments[i]]]
    Cc <- average(Ce[Seg_i[[length(Seg_i)]], , drop = FALSE])
    LinearDist[i] <- sqrt(sum((Cc - BranchBase)^2))
  }
  
  ord <- order(LinearDist, decreasing = TRUE)
  LinearDist <- LinearDist[ord]
  Segments   <- Segments[ord]
  
  # Definisce i rami candidati (fusioni successive)
  Branches <- vector("list", t)
  
  for (j in seq_len(t)) {
    Seg  <- Segs[[BI]]
    spar <- SPar
    
    if (Segments[j] != BI) {
      SubSegs <- integer(200)
      SubSegs[1] <- Segments[j]
      k <- 1L
      S <- Segments[j]
      safety_chain <- 0L
      while (S != BI) {
        S <- SPar[S, 1]
        k <- k + 1L
        SubSegs[k] <- S
        safety_chain <- safety_chain + 1L
        if (safety_chain > ns + 10L) {
          warning("search_branch_top: risalita verso BI non converge (probabile ciclo in SPar) - uso fallback BranchTop=BI")
          return(BI)
        }
      }
      
      Segment <- vector("list", BUFLEN)
      a <- length(Seg)
      Segment[seq_len(a)] <- Seg
      a <- a + 1L
      
      if (k > 2L) {
        for (i2 in seq_len(k - 2L)) {
          I_seg <- SubSegs[k - i2]
          J_seg <- SubSegs[k - i2 - 1L]
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
          SubSegs[k - i2] <- 1L
        }
      }
      
      I_seg <- SubSegs[1]
      SP <- spar[I_seg, 2]
      SegC <- Segs[[I_seg]]
      L_len <- length(SegC)
      if (SP >= a - 2L) {
        Segment[a:(a + L_len - 1L)] <- SegC
        a <- a + L_len - 1L
      } else {
        Segment[(SP + 1L):(SP + L_len)] <- SegC
        a <- SP + L_len
      }
      Branches[[j]] <- Segment[seq_len(a)]
    } else {
      Branches[[j]] <- Seg
    }
  }
  
  # Calcola le lunghezze approssimate dei rami candidati, in ordine di
  # distanza lineare decrescente, fermandosi al primo che risulta
  # abbastanza dritto
  N <- ceiling(0.25 / dmin / 1.4)
  i <- 1L
  Continue <- TRUE
  Lengths <- numeric(t)
  BranchTop <- Segments[1]  # fallback difensivo, non nell'originale
  
  while (i <= t && Continue) {
    Seg2 <- Branches[[i]]
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
      Lengths[i] <- sum(sqrt(rowSums(V * V)))
    } else {
      Lengths[i] <- 0
    }
    
    ratio_i <- if (LinearDist[i] > 0) Lengths[i] / LinearDist[i] else Inf
    if (is.finite(ratio_i) && ratio_i < 1.20 && LinearDist[i] > 0.75 * LinearDist[1]) {
      Continue <- FALSE
      BranchTop <- Segments[i]
    }
    i <- i + 1L
  }
  
  # Se nessun candidato ha soddisfatto il criterio stretto, riprova con
  # una soglia piu' permissiva; in ultima istanza usa quello con la
  # distanza lineare maggiore
  if (Continue) {
    Lr <- Lengths / LinearDist
    Lr[!is.finite(Lr)] <- Inf  # NaN (0/0) trattato come "molto curvo", non soddisfa mai il criterio < 1.4
    i <- 1L
    while (i <= t && Lr[i] > 1.4 && LinearDist[i] > 0.75 * LinearDist[1]) {
      i <- i + 1L
    }
    if (i <= t) {
      BranchTop <- Segments[i]
    } else {
      BranchTop <- Segments[1]
    }
  }
  
  BranchTop
}