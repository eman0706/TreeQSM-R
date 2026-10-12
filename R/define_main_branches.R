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

#' Pairwise Distances (pdist2)
#'
#' Equivalente di pdist2() di MATLAB - matrice delle distanze euclidee
#' tra le righe di A e le righe di B.
#' @param A Matrice n x d
#' @param B Matrice m x d
#' @return Matrice n x m di distanze
#' @export
pdist2 <- function(A, B) {
  A <- matrix(A, ncol = ncol(as.matrix(B)))
  B <- as.matrix(B)
  sq <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  sqrt(pmax(sq, 0))
}

#' Define Main Branches
#'
#' Porting fedele di define_main_branches (sotto-funzione di
#' tree_sets.m, Pasi Raumonen). Se esiste una segmentazione precedente
#' (dal primo passaggio), la usa per rendere connessi i cover set nei
#' rami principali (tronco e rami di ordine 1-3) - ripara la
#' connettivita' tramite ricerca spaziale cubica (cubical_partition),
#' unendo esplicitamente le componenti sconnesse al vicino piu' vicino
#' della stessa categoria di ramo principale.
#'
#' NOTA DI FEDELTA': l'ultima chiamata a connected_components()
#' nell'originale usa solo 3 argomenti (Nei, Trunk, Fal), a differenza
#' delle altre nella stessa funzione (4 argomenti, con MinSize=1
#' esplicito) - senza il vero sorgente di connected_components.m per
#' verificare la gestione di argomenti omessi, si assume MinSize=1
#' anche qui (coerente con lo scopo evidente: classificare tutte le
#' componenti, nessuna scartata per dimensione minima).
#'
#' NOTA DI SICUREZZA: aggiunto un limite alle iterazioni del ciclo
#' "while (n_comp > 1)" di riparazione connettivita' (non presente
#' nell'originale) - lezione della sessione: un ciclo di
#' riparazione/fusione iterativa puo' non convergere mai su dati
#' patologici, come visto con modify_topology.
#'
#' @param cover Cover set (ball, neighbor, center)
#' @param segment Segmentazione precedente (deve avere SegmentOfPoint,
#'   branch1indexes, branch2indexes, branch3indexes - da
#'   correct_segments())
#' @param aux Lista ausiliaria (Ce, Fal, Ind, nb)
#' @param inputs Parametri (usa PatchDiam2Max)
#' @return Lista con Trunk (vettore logico) e cover (con neighbor aggiornato)
#' @export
define_main_branches <- function(cover, segment, aux, inputs) {
  Bal <- cover$ball
  Nei <- cover$neighbor
  Ce  <- aux$Ce
  nb  <- length(Bal)
  
  MainBranches <- integer(nb)
  SegmentOfPoint <- segment$SegmentOfPoint
  max_seg <- suppressWarnings(max(SegmentOfPoint, na.rm = TRUE))
  if (!is.finite(max_seg) || max_seg < 1) {
    return(list(Trunk = aux$Fal, cover = cover))
  }
  
  MainBranchIndexes <- rep(FALSE, max_seg)
  MainBranchIndexes[1] <- TRUE
  b1 <- segment$branch1indexes; b1 <- b1[!is.na(b1) & b1 >= 1 & b1 <= max_seg]
  b2 <- segment$branch2indexes; b2 <- b2[!is.na(b2) & b2 >= 1 & b2 <= max_seg]
  b3 <- segment$branch3indexes; b3 <- b3[!is.na(b3) & b3 >= 1 & b3 <= max_seg]
  if (length(b1) > 0) MainBranchIndexes[b1] <- TRUE
  if (length(b2) > 0) MainBranchIndexes[b2] <- TRUE
  if (length(b3) > 0) MainBranchIndexes[b3] <- TRUE
  
  for (i in seq_len(nb)) {
    pts <- Bal[[i]]
    if (length(pts) == 0) next
    BranchInd <- SegmentOfPoint[pts]
    BranchInd <- BranchInd[!is.na(BranchInd) & BranchInd != 0]
    if (length(BranchInd) > 0) {
      ind <- min(BranchInd)
      if (ind >= 1 && ind <= max_seg && MainBranchIndexes[ind]) {
        MainBranches[i] <- min(BranchInd)
      }
    }
  }
  
  Trunk <- aux$Fal
  Trunk[MainBranches > 0] <- TRUE
  
  ## Aggiorna i vicini per rendere connessi i rami principali
  cp <- cubical_partition(Ce, 3 * inputs$PatchDiam2Max, 10)
  Par <- cp$Partition
  CC  <- cp$CubeCoord
  BI  <- max(MainBranches)
  N_dim <- dim(Par)
  
  if (BI >= 1) {
    for (i in seq_len(BI)) {
      if (i > length(MainBranchIndexes) || !MainBranchIndexes[i]) next
      Branch <- MainBranches == i
      cc_res <- connected_components(Nei, Branch, 1, aux$Fal)
      Comps <- cc_res$Components
      n_comp <- length(Comps)
      
      safety <- 0L
      safety_limit <- nb + 100L
      while (n_comp > 1) {
        safety <- safety + 1L
        if (safety > safety_limit) {
          warning("define_main_branches: riparazione connettivita' non converge - interrompo con lo stato corrente")
          break
        }
        for (j in seq_len(n_comp)) {
          comp <- Comps[[j]]
          NC <- length(comp)
          c_coords <- unique(CC[comp, , drop = FALSE])
          m <- nrow(c_coords)
          t_val <- 0L
          NearSets <- integer(0)
          
          safety2 <- 0L
          while (length(NearSets) == 0) {
            safety2 <- safety2 + 1L
            if (safety2 > max(N_dim) + 5L) break  # esaurita la griglia, nessun vicino trovabile
            NearSets_bool <- aux$Fal
            t_val <- t_val + 1L
            for (k in seq_len(m)) {
              x1 <- max(1, c_coords[k, 1] - t_val); x2 <- min(c_coords[k, 1] + t_val, N_dim[1])
              y1 <- max(1, c_coords[k, 2] - t_val); y2 <- min(c_coords[k, 2] + t_val, N_dim[2])
              z1 <- max(1, c_coords[k, 3] - t_val); z2 <- min(c_coords[k, 3] + t_val, N_dim[3])
              balls0 <- Par[x1:x2, y1:y2, z1:z2]
              balls <- unlist(balls0, use.names = FALSE)
              if (length(balls) > 0) {
                balls <- balls[balls >= 1 & balls <= nb]
                if (length(balls) > 0) {
                  bb <- balls[Branch[balls]]
                  if (length(bb) > 0) NearSets_bool[bb] <- TRUE
                }
              }
            }
            NearSets_bool[comp] <- FALSE
            NearSets <- aux$Ind[NearSets_bool]
          }
          
          if (length(NearSets) > 0) {
            d_mat <- pdist2(Ce[comp, , drop = FALSE], Ce[NearSets, , drop = FALSE])
            if (NC == 1 && length(NearSets) == 1) {
              IU <- 1L; JU <- 1L
            } else if (NC == 1) {
              JU <- which.min(d_mat[1, ]); IU <- 1L
            } else if (length(NearSets) == 1) {
              IU <- which.min(d_mat[, 1]); JU <- 1L
            } else {
              col_mins <- apply(d_mat, 2, min)
              IU_per_col <- apply(d_mat, 2, which.min)
              JU <- which.min(col_mins)
              IU <- IU_per_col[JU]
            }
            I_conn <- comp[IU]
            J_conn <- NearSets[JU]
            Nei[[I_conn]] <- c(Nei[[I_conn]], J_conn)
            Nei[[J_conn]] <- c(Nei[[J_conn]], I_conn)
          }
        }
        cc_res <- connected_components(Nei, Branch, 1, aux$Fal)
        Comps <- cc_res$Components
        n_comp <- length(Comps)
      }
    }
  }
  
  ## Collega i rami di primo ordine al tronco
  Stem <- aux$Ind[MainBranches == 1]
  MainBranchIndexes2 <- rep(FALSE, max_seg)
  b1b <- segment$branch1indexes; b1b <- b1b[!is.na(b1b) & b1b >= 1 & b1b <= max_seg]
  if (length(b1b) > 0) MainBranchIndexes2[b1b] <- TRUE
  BI2 <- if (length(segment$branch1indexes) > 0) {
    suppressWarnings(max(segment$branch1indexes, na.rm = TRUE))
  } else {
    0L
  }
  if (!is.finite(BI2)) BI2 <- 0L
  
  if (BI2 >= 2) {
    for (i in 2:BI2) {
      if (i > length(MainBranchIndexes2) || !MainBranchIndexes2[i]) next
      BranchSet <- aux$Ind[MainBranches == i]
      if (length(BranchSet) == 0) next
      
      neigh_all <- unlist(Nei[BranchSet], use.names = FALSE)
      neigh_all <- neigh_all[!is.na(neigh_all) & neigh_all >= 1 & neigh_all <= nb]
      Neighbors_is1 <- if (length(neigh_all) > 0) MainBranches[neigh_all] == 1 else logical(0)
      
      if (!any(Neighbors_is1)) {
        d_mat <- pdist2(Ce[BranchSet, , drop = FALSE], Ce[Stem, , drop = FALSE])
        if (length(BranchSet) > 1 && length(Stem) > 1) {
          col_mins <- apply(d_mat, 2, min)
          IU_per_col <- apply(d_mat, 2, which.min)
          J_idx <- which.min(col_mins)
          I_idx <- IU_per_col[J_idx]
        } else if (length(BranchSet) == 1 && length(Stem) > 1) {
          J_idx <- which.min(d_mat[1, ]); I_idx <- 1L
        } else if (length(Stem) == 1 && length(BranchSet) > 1) {
          I_idx <- which.min(d_mat[, 1]); J_idx <- 1L
        } else if (length(BranchSet) == 1 && length(Stem) == 1) {
          I_idx <- 1L; J_idx <- 1L
        } else {
          I_idx <- NA; J_idx <- NA
        }
        if (!is.na(I_idx) && !is.na(J_idx)) {
          I_conn <- BranchSet[I_idx]
          J_conn <- Stem[J_idx]
          Nei[[I_conn]] <- c(Nei[[I_conn]], J_conn)
          Nei[[J_conn]] <- c(Nei[[J_conn]], I_conn)
        }
      }
    }
  }
  cover$neighbor <- Nei
  
  ## Verifica che il tronco resti in un'unica componente; se no, sceglie
  ## quella che contiene il vero stem (branch1)
  cc_final <- connected_components(cover$neighbor, Trunk, 1, aux$Fal)
  comps_final <- cc_final$Components
  cs_final <- cc_final$CompSize
  if (length(cs_final) > 1) {
    ord <- order(cs_final, decreasing = TRUE)
    cs_final <- cs_final[ord]
    comps_final <- comps_final[ord]
    Stem_bool <- MainBranches == 1
    Trunk <- aux$Fal
    ii <- 1L
    C_comp <- comps_final[[ii]]
    while (ii <= length(cs_final) && !any(Stem_bool[C_comp])) {
      ii <- ii + 1L
      if (ii > length(cs_final)) break
      C_comp <- comps_final[[ii]]
    }
    if (ii <= length(cs_final)) Trunk[C_comp] <- TRUE
  }
  
  list(Trunk = Trunk, cover = cover)
}
