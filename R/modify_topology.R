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

#' Modify Topology
#'
#' Porting fedele di modify_topology (sotto-funzione di
#' correct_segments.m, Pasi Raumonen). Rende il tronco e i rami "il piu'
#' lunghi possibile", fondendo segmenti che sono in realta' la
#' continuazione fisica dello stesso ramo (invece di restare frammentati
#' come li produce segments.m, che deliberatamente taglia troppo,
#' sapendo che questa funzione ripulira' subito dopo). Processa i
#' segmenti in ordine crescente di ordine di ramificazione (prima il
#' tronco, poi tutti i rami di 1' ordine, poi 2', ecc.), cercando per
#' ciascuno la vera punta (tramite search_stem_top/search_branch_top
#' per tronco/1' ordine, o per distanza lineare diretta per ordini
#' superiori) e fondendo la catena di segmenti intermedi.
#'
#' NOTA DI IMPLEMENTAZIONE: SPar e' mantenuta come matrice (per
#' compatibilita' con il resto della pipeline) pre-allocata
#' generosamente (2x il numero iniziale di segmenti, stessa euristica
#' del buffer "Fal = false(2*ns,1)" dell'originale) invece di essere
#' ridimensionata ad ogni nuovo segmento creato durante le divisioni.
#' Segs e SChi restano liste R, che crescono naturalmente.
#'
#' NOTA FEDELTA' (trovata durante l'autoverifica): le soglie di
#' confronto "SP >= N-2" (nel ciclo di fusione intermedio) e
#' "SP >= N-3" (nella combinazione finale dell'ultimo segmento) sono
#' DIVERSE tra loro nell'originale - non un refuso, preservate cosi'
#' come sono.
#'
#' @param P Point cloud
#' @param Ce Centri dei cover set
#' @param Bal cover$ball
#' @param Segs segment$segments
#' @param SPar segment$ParentSegment (matrice ns x 2)
#' @param SChi segment$ChildSegment
#' @param dmin Diametro minimo di riferimento
#' @return Lista con Segs, SPar, SChi aggiornati
#' @export
modify_topology <- function(P, Ce, Bal, Segs, SPar, SChi, dmin) {
  ns0 <- length(Segs)
  
  # SPar pre-allocata generosamente (2x, come "Fal = false(2*ns,1)"
  # nell'originale suggerisce), troncata alla fine.
  SPar_buf <- matrix(0L, 2L * ns0, 2L)
  SPar_buf[seq_len(ns0), ] <- SPar[seq_len(ns0), ]
  SPar <- SPar_buf
  
  ns <- ns0
  SubSegments <- integer(max(200L, ceiling(ns0 / 5) + 10L))
  SegInd <- 1L
  UnMod <- rep(TRUE, ns0)
  UnMod[SegInd] <- FALSE
  BranchOrder <- 0L
  ChildSegInd <- 1L
  SegChildren <- integer(0)
  safety_main <- 0L
  safety_main_limit <- ns0 * 20L + 1000L
  
  while (any(UnMod[seq_len(ns)])) {
    safety_main <- safety_main + 1L
    if (safety_main > safety_main_limit) {
      warning("modify_topology: ciclo principale non converge dopo troppe iterazioni - interrompo con lo stato corrente")
      break
    }
    ChildSegs <- SChi[[SegInd]]
    
    if (length(Segs[[SegInd]]) > 0 && length(ChildSegs) > 0) {
      
      if (SegInd > 1L && BranchOrder > 1L) {
        # Rami di 2' ordine e superiori: cerca la punta del sotto-ramo
        # con maggiore distanza lineare dalla base del ramo corrente
        SubSegments[1] <- SegInd
        NSubSegs <- 2L
        ChSegs2 <- ChildSegs
        n_last <- 0L
        while (length(ChSegs2) > 0) {
          n_last <- length(ChSegs2)
          SubSegments[NSubSegs:(NSubSegs + n_last - 1L)] <- ChSegs2
          ChSegs2 <- unlist(SChi[ChSegs2])
          NSubSegs <- NSubSegs + n_last
        }
        if (NSubSegs > 2L) {
          NSubSegs <- NSubSegs - n_last
        }
        
        Top <- matrix(0, NSubSegs, 3)
        for (i in seq_len(NSubSegs)) {
          seg_i <- Segs[[SubSegments[i]]]
          Top[i, ] <- Ce[seg_i[[length(seg_i)]][1], ]
        }
        
        BotLayer <- Segs[[SegInd]][[1]]
        Bottom <- average(Ce[BotLayer, , drop = FALSE])
        
        V <- Top - matrix(Bottom, NSubSegs, 3, byrow = TRUE)
        d <- rowSums(V * V)
        I_max <- which.max(d)
        TipSeg <- SubSegments[I_max[1]]
        
      } else if (SegInd > 1L && BranchOrder <= 1L) {
        # Rami di 1' ordine
        TipSeg <- search_branch_top(P, Ce, Bal, Segs, SPar, SChi, dmin, SegInd)
      } else {
        # Tronco
        TipSeg <- search_stem_top(P, Ce, Bal, Segs, SPar, dmin)
      }
      
      if (TipSeg != SegInd) {
        SubSegments[1] <- TipSeg
        NSubSegs <- 1L
        safety_chain <- 0L
        while (TipSeg != SegInd) {
          TipSeg <- SPar[TipSeg, 1]
          NSubSegs <- NSubSegs + 1L
          SubSegments[NSubSegs] <- TipSeg
          safety_chain <- safety_chain + 1L
          if (safety_chain > ns0 * 2L + 10L) {
            warning("modify_topology: risalita verso SegInd non converge (probabile ciclo in SPar) - salto la fusione per questo segmento")
            TipSeg <- SegInd
            break
          }
        }
        
        if (NSubSegs > 2L) {
          for (i in seq_len(NSubSegs - 2L)) {
            I_seg <- SubSegments[NSubSegs - i]
            J_seg <- SubSegments[NSubSegs - i - 1L]
            SP <- SPar[I_seg, 2]
            SegP <- Segs[[SegInd]]
            SegC <- Segs[[I_seg]]
            N_len <- length(SegP)
            sp <- SPar[J_seg, 2]
            
            if (SP >= N_len - 2L) {
              # Usa l'intero genitore
              Segs[[SegInd]] <- c(SegP, SegC[seq_len(sp)])
              if (sp < length(SegC)) {
                # usa solo parte del segmento figlio
                Segs[[I_seg]] <- SegC[(sp + 1L):length(SegC)]
                SPar[I_seg, 2] <- N_len + sp
                
                ChildSegs2 <- SChi[[I_seg]]
                K <- SPar[ChildSegs2, 2] <= sp
                c_keep <- ChildSegs2[!K]
                SChi[[I_seg]] <- c_keep
                if (length(c_keep) > 0) SPar[c_keep, 2] <- SPar[c_keep, 2] - sp
                ChildSegs2 <- ChildSegs2[K]
                SChi[[SegInd]] <- c(SChi[[SegInd]], ChildSegs2)
                if (length(ChildSegs2) > 0) {
                  SPar[ChildSegs2, 1] <- SegInd
                  SPar[ChildSegs2, 2] <- N_len + SPar[ChildSegs2, 2]
                }
              } else {
                # usa l'intero segmento figlio
                Segs[[I_seg]] <- list()
                SPar[I_seg, 1] <- 0L
                UnMod[I_seg] <- FALSE
                
                ChildSegs2 <- SChi[[I_seg]]
                SChi[[I_seg]] <- integer(0)
                c_diff <- setdiff(SChi[[SegInd]], I_seg)
                SChi[[SegInd]] <- c(c_diff, ChildSegs2)
                if (length(ChildSegs2) > 0) {
                  SPar[ChildSegs2, 1] <- SegInd
                  SPar[ChildSegs2, 2] <- N_len + SPar[ChildSegs2, 2]
                }
              }
              SubSegments[NSubSegs - i] <- SegInd
              
            } else {
              # Divide il segmento genitore in due parti
              ns <- ns + 1L
              if (ns > nrow(SPar)) SPar <- rbind(SPar, matrix(0L, nrow(SPar), 2L))
              Segs[[ns]] <- SegP[(SP + 1L):length(SegP)]
              SPar[ns, 1] <- SegInd
              SPar[ns, 2] <- SP
              UnMod[ns] <- TRUE
              
              Segs[[SegInd]] <- c(SegP[seq_len(SP)], SegC[seq_len(sp)])
              
              ChildSegs2 <- SChi[[SegInd]]
              K <- SPar[ChildSegs2, 2] > SP
              SChi[[SegInd]] <- ChildSegs2[!K]
              ChildSegs2 <- ChildSegs2[K]
              SChi[[ns]] <- ChildSegs2
              if (length(ChildSegs2) > 0) {
                SPar[ChildSegs2, 1] <- ns
                SPar[ChildSegs2, 2] <- SPar[ChildSegs2, 2] - SP
              }
              SChi[[SegInd]] <- c(SChi[[SegInd]], ns)
              
              if (sp < length(SegC)) {
                Segs[[I_seg]] <- SegC[(sp + 1L):length(SegC)]
                SPar[I_seg, 2] <- SP + sp
                
                ChildSegs2 <- SChi[[I_seg]]
                K <- SPar[ChildSegs2, 2] <= sp
                SChi[[I_seg]] <- ChildSegs2[!K]
                if (any(!K)) SPar[ChildSegs2[!K], 2] <- SPar[ChildSegs2[!K], 2] - sp
                ChildSegs2 <- ChildSegs2[K]
                SChi[[SegInd]] <- c(SChi[[SegInd]], ChildSegs2)
                if (length(ChildSegs2) > 0) {
                  SPar[ChildSegs2, 1] <- SegInd
                  SPar[ChildSegs2, 2] <- SP + SPar[ChildSegs2, 2]
                }
              } else {
                Segs[[I_seg]] <- list()
                SPar[I_seg, 1] <- 0L
                UnMod[I_seg] <- FALSE
                
                ChildSegs2 <- SChi[[I_seg]]
                c_diff <- setdiff(SChi[[SegInd]], I_seg)
                SChi[[SegInd]] <- c(c_diff, ChildSegs2)
                if (length(ChildSegs2) > 0) {
                  SPar[ChildSegs2, 1] <- SegInd
                  SPar[ChildSegs2, 2] <- SP + SPar[ChildSegs2, 2]
                }
              }
              SubSegments[NSubSegs - i] <- SegInd
            }
          }
        }
        
        # Combina l'ultimo segmento della catena al ramo
        # NOTA FEDELTA': qui la soglia e' N-3, diversa da N-2 usata sopra
        I_seg <- SubSegments[1]
        SP <- SPar[I_seg, 2]
        SegP <- Segs[[SegInd]]
        SegC <- Segs[[I_seg]]
        N_len <- length(SegP)
        
        if (SP >= N_len - 3L) {
          Segs[[SegInd]] <- c(SegP, SegC)
          Segs[[I_seg]] <- list()
          SPar[I_seg, 1] <- 0L
          UnMod[I_seg] <- FALSE
          
          ChildSegs2 <- SChi[[I_seg]]
          c_diff <- setdiff(SChi[[SegInd]], I_seg)
          SChi[[SegInd]] <- c(c_diff, ChildSegs2)
          if (length(ChildSegs2) > 0) {
            SPar[ChildSegs2, 1] <- SegInd
            SPar[ChildSegs2, 2] <- N_len + SPar[ChildSegs2, 2]
          }
        } else {
          ns <- ns + 1L
          if (ns > nrow(SPar)) SPar <- rbind(SPar, matrix(0L, nrow(SPar), 2L))
          Segs[[ns]] <- SegP[(SP + 1L):length(SegP)]
          SPar[ns, ] <- c(SegInd, SP)
          Segs[[SegInd]] <- c(SegP[seq_len(SP)], SegC)
          Segs[[I_seg]] <- list()
          SPar[I_seg, 1] <- 0L
          UnMod[ns] <- TRUE
          UnMod[I_seg] <- FALSE
          
          ChildSegs2 <- SChi[[SegInd]]
          K <- SPar[ChildSegs2, 2] > SP
          # FIX: l'originale aggiunge "ns" qui stesso (SChi{SegInd} =
          # [ChildSegs(~K); ns]) - mancava in questa traduzione, lasciando
          # "ns" orfano: SPar[ns,1] lo dichiarava figlio di SegInd, ma
          # SChi[[SegInd]] non lo conteneva mai, rendendolo irraggiungibile
          # dalla scansione top-down di modify_topology (causa del loop
          # che non convergeva mai, trovato con la diagnostica mirata).
          SChi[[SegInd]] <- c(ChildSegs2[!K], ns)
          ChildSegs2b <- ChildSegs2[K]
          SChi[[ns]] <- ChildSegs2b
          if (length(ChildSegs2b) > 0) {
            SPar[ChildSegs2b, 1] <- ns
            SPar[ChildSegs2b, 2] <- SPar[ChildSegs2b, 2] - SP
          }
          
          ChildSegs2 <- SChi[[I_seg]]
          c_diff <- setdiff(SChi[[SegInd]], I_seg)
          SChi[[SegInd]] <- c(c_diff, ChildSegs2)
          if (length(ChildSegs2) > 0) {
            SPar[ChildSegs2, 1] <- SegInd
            SPar[ChildSegs2, 2] <- SP + SPar[ChildSegs2, 2]
          }
        }
      }
      UnMod[SegInd] <- FALSE
    } else {
      UnMod[SegInd] <- FALSE
    }
    
    # Seleziona il prossimo segmento, seguendo l'ordine crescente di
    # ramificazione
    if (BranchOrder > 0L && length(SegChildren) > 0 && any(UnMod[SegChildren])) {
      ChildSegInd <- ChildSegInd + 1L
      SegInd <- SegChildren[ChildSegInd]
    } else if (BranchOrder == 0L) {
      BranchOrder <- BranchOrder + 1L
      SegChildren <- SChi[[1]]
      if (length(SegChildren) > 0) {
        SegInd <- SegChildren[1]
      } else {
        UnMod[seq_len(ns)] <- FALSE
      }
    } else {
      BranchOrder <- BranchOrder + 1L
      i <- 1L
      SegChildren <- SChi[[1]]
      while (i < BranchOrder && length(SegChildren) > 0) {
        i <- i + 1L
        L_len <- sapply(SChi[SegChildren], length)
        Keep <- L_len > 0
        SegChildren <- SegChildren[Keep]
        SegChildren <- unlist(SChi[SegChildren])
      }
      if (length(SegChildren) > 0) {
        I_unmod <- UnMod[SegChildren]
        if (any(I_unmod)) {
          SegChildren <- SegChildren[I_unmod]
          SegInd <- SegChildren[1]
          ChildSegInd <- 1L
        }
      }
    }
  }
  
  # Rimuove i segmenti vuoti e ricompatta gli indici
  Empty <- rep(TRUE, ns)
  for (i in seq_len(ns)) {
    if (length(Segs[[i]]) == 0) Empty[i] <- FALSE
  }
  Segs <- Segs[Empty]
  Ind <- seq_len(ns)
  n_keep <- sum(Empty)
  Ind[Empty] <- seq_len(n_keep)
  SPar <- SPar[seq_len(ns), , drop = FALSE][Empty, , drop = FALSE]
  J <- SPar[, 1] > 0
  if (any(J)) SPar[J, 1] <- Ind[SPar[J, 1]]
  
  SChi_kept <- SChi[seq_len(ns)]
  for (i in seq_len(ns)) {
    if (Empty[i]) {
      ChildSegs2 <- SChi_kept[[i]]
      if (length(ChildSegs2) > 0) {
        SChi_kept[[i]] <- Ind[ChildSegs2]
      }
    }
  }
  SChi <- SChi_kept[Empty]
  ns <- n_keep
  
  list(Segs = Segs, SPar = SPar, SChi = SChi)
}