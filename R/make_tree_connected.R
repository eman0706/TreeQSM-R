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

#' Make Tree Connected
#'
#' Porting fedele di make_tree_connected (sotto-funzione di
#' tree_sets.m, Pasi Raumonen). Aggiorna la relazione di vicinato in
#' modo che l'intero albero sia un'unica componente connessa. Espande
#' il Trunk il piu' possibile, poi esamina i cover set rimasti
#' isolati ("Other") componente per componente, classificando ciascuno
#' come:
#'   - rumore da escludere (Forb), se troppo lontano da tutto, o
#'     piccolo e a distanza moderata dal Trunk;
#'   - parte del Trunk, se abbastanza vicino;
#'   - lasciato isolato, se la distanza supera la finestra di ricerca
#'     corrente (che si allarga ad ogni round finche' tutto e' classificato).
#'
#' Questo e' il meccanismo geometrico (non un filtro statistico sui
#' punti grezzi) che nella sessione si sospetta possa correttamente
#' escludere contaminazione da vegetazione bassa/sottobosco erroneamente
#' inclusa nella nuvola di punti durante l'estrazione LiDAR - il
#' problema diagnosticato empiricamente su piu' alberi.
#'
#' NOTA DI SICUREZZA: aggiunti limiti alle iterazioni dei cicli "while"
#' (non presenti nell'originale) - lezione della sessione: cicli di
#' classificazione/espansione iterativa possono non convergere mai su
#' dati patologici.
#'
#' @param cover Cover set (neighbor, center - via aux$Ce)
#' @param aux Lista ausiliaria (Ce, Fal, Ind, nb)
#' @param Forb Vettore logico dei cover set proibiti/esclusi
#' @param Base Indici dei cover set della base
#' @param Trunk Vettore logico dei cover set del tronco/rami principali
#' @param inputs Parametri (usa PatchDiam1, OnlyTree)
#' @return Lista con cover (neighbor aggiornato) e Forb aggiornato
#' @export
make_tree_connected <- function(cover, aux, Forb, Base, Trunk, inputs) {
  Nei <- cover$neighbor
  Ce  <- aux$Ce
  nb  <- aux$nb
  
  ## Espande il tronco il piu' possibile
  Trunk[Forb] <- FALSE
  Exp <- Trunk
  safety_expand <- 0L
  while (any(Exp)) {
    safety_expand <- safety_expand + 1L
    if (safety_expand > nb + 100L) {
      warning("make_tree_connected: espansione del tronco non converge - interrompo")
      break
    }
    exp_idx <- which(Exp)
    new_exp <- unlist(Nei[exp_idx], use.names = FALSE)
    new_exp <- new_exp[!is.na(new_exp) & new_exp >= 1 & new_exp <= nb]
    Exp2 <- aux$Fal
    if (length(new_exp) > 0) Exp2[new_exp] <- TRUE
    Exp2[Trunk] <- FALSE
    Exp2[Forb] <- FALSE
    Exp2[Base] <- FALSE
    Trunk[Exp2] <- TRUE
    Exp <- Exp2
  }
  
  ## Definisce "Other": cover set non ancora connessi a Trunk o Forb
  Other <- !aux$Fal
  Other[Forb] <- FALSE
  Other[Trunk] <- FALSE
  Other[Base] <- FALSE
  
  ## Parametri dell'estensione dello spazio di ricerca e della
  ## dimensione minima accettabile di componente
  k0 <- min(10, ceiling(0.2 / inputs$PatchDiam1))
  k <- k0
  if (isTRUE(inputs$OnlyTree)) {
    Cmin <- 0
  } else {
    Cmin <- ceiling(0.1 / inputs$PatchDiam1)
  }
  
  if (any(Other)) {
    cc0 <- connected_components(Nei, Other, 1, aux$Fal)
    Comps <- cc0$Components
    nc <- length(Comps)
    NonClassified <- rep(TRUE, nc)
  } else {
    NonClassified <- FALSE
    nc <- 0L
    Comps <- list()
  }
  
  bottom <- min(Ce[Base, 3])
  
  safety_outer <- 0L
  safety_outer_limit <- 1000L
  while (any(NonClassified)) {
    safety_outer <- safety_outer + 1L
    if (safety_outer > safety_outer_limit) {
      warning("make_tree_connected: ciclo esterno non converge - interrompo con lo stato corrente")
      break
    }
    npre <- sum(NonClassified)
    again <- TRUE
    
    cp <- cubical_partition(Ce, k * inputs$PatchDiam1)
    Par <- cp$Partition
    CC  <- cp$CubeCoord
    N_dim <- dim(Par)
    Neighbors <- vector("list", nc)
    Sizes <- matrix(0L, nc, 2)
    Pass <- rep(TRUE, nc)
    first_round <- TRUE
    
    safety_inner <- 0L
    safety_inner_limit <- 1000L
    while (again) {
      safety_inner <- safety_inner + 1L
      if (safety_inner > safety_inner_limit) {
        warning("make_tree_connected: ciclo interno non converge - interrompo")
        break
      }
      for (i in seq_len(nc)) {
        if (!NonClassified[i] || !Pass[i]) next
        comp <- Comps[[i]]
        
        neigh_comp <- unlist(Nei[comp], use.names = FALSE)
        neigh_comp <- neigh_comp[!is.na(neigh_comp) & neigh_comp >= 1 & neigh_comp <= nb]
        J_forb <- if (length(neigh_comp) > 0) Forb[neigh_comp] else logical(0)
        
        if (any(J_forb)) {
          NonClassified[i] <- FALSE
          Forb[comp] <- TRUE
          Other[comp] <- FALSE
        } else {
          NC <- length(comp)
          if (first_round) {
            c_coords <- unique(CC[comp, , drop = FALSE])
            m <- nrow(c_coords)
            B_list <- vector("list", m)
            for (j in seq_len(m)) {
              xr <- (c_coords[j, 1] - 1):(c_coords[j, 1] + 1)
              yr <- (c_coords[j, 2] - 1):(c_coords[j, 2] + 1)
              zr <- (c_coords[j, 3] - 1):(c_coords[j, 3] + 1)
              xr <- xr[xr >= 1 & xr <= N_dim[1]]
              yr <- yr[yr >= 1 & yr <= N_dim[2]]
              zr <- zr[zr >= 1 & zr <= N_dim[3]]
              if (length(xr) == 0 || length(yr) == 0 || length(zr) == 0) {
                B_list[[j]] <- integer(0)
              } else {
                balls <- Par[xr, yr, zr]
                B_list[[j]] <- unlist(balls, use.names = FALSE)
              }
            }
            NearSets <- unlist(B_list, use.names = FALSE)
            in_comp <- aux$Fal
            in_comp[comp] <- TRUE
            if (length(NearSets) > 0) {
              NearSets <- NearSets[!in_comp[NearSets]]
              NearSets <- unique(NearSets)
            }
            Neighbors[[i]] <- NearSets
            if (length(NearSets) == 0) {
              Pass[i] <- FALSE
            }
            if (length(NearSets) > 0) NearSets <- NearSets[!Other[NearSets]]
          } else {
            NearSets <- Neighbors[[i]]
            if (length(NearSets) > 0) NearSets <- NearSets[!Other[NearSets]]
          }
          
          I_trunk <- if (length(NearSets) > 0) Trunk[NearSets] else logical(0)
          J_forb2 <- if (length(NearSets) > 0) Forb[NearSets] else logical(0)
          trunk_set <- NearSets[I_trunk]
          forb_set  <- NearSets[J_forb2]
          
          if (length(trunk_set) != Sizes[i, 1] || length(forb_set) != Sizes[i, 2]) {
            Sizes[i, ] <- c(length(trunk_set), length(forb_set))
            
            if (NC > 100) {
              hmin <- min(Ce[comp, 3])
              H <- max(Ce[comp, 3]) - hmin
              if (H > 5 && hmin < bottom + 5) {
                I_h <- Ce[NearSets, 3] < hmin + 0.5
                NearSets <- NearSets[I_h]
                I_trunk <- if (length(NearSets) > 0) Trunk[NearSets] else logical(0)
                J_forb2 <- if (length(NearSets) > 0) Forb[NearSets] else logical(0)
                trunk_set <- NearSets[I_trunk]
                forb_set  <- NearSets[J_forb2]
              }
            }
            
            if (length(trunk_set) > 0) {
              d_mat <- pdist2(Ce[comp, , drop = FALSE], Ce[trunk_set, , drop = FALSE])
              if (NC == 1 && length(trunk_set) == 1) {
                dt <- d_mat[1, 1]; IC <- 1L; IT <- 1L
              } else if (NC == 1) {
                IT <- which.min(d_mat[1, ]); dt <- d_mat[1, IT]; IC <- 1L
              } else if (length(trunk_set) == 1) {
                IC <- which.min(d_mat[, 1]); dt <- d_mat[IC, 1]; IT <- 1L
              } else {
                col_mins <- apply(d_mat, 2, min)
                IC_per_col <- apply(d_mat, 2, which.min)
                IT <- which.min(col_mins)
                dt <- col_mins[IT]
                IC <- IC_per_col[IT]
              }
            } else {
              dt <- 700
              IC <- NA_integer_; IT <- NA_integer_
            }
            
            if (length(forb_set) > 0) {
              d_mat_f <- pdist2(Ce[comp, , drop = FALSE], Ce[forb_set, , drop = FALSE])
              df <- min(d_mat_f)
            } else {
              df <- 1000
            }
            
            if ((dt > 12 && dt < 100) || (NC < Cmin && dt > 0.5 && dt < 10)) {
              # Rimuove il piccolo componente isolato
              Forb[comp] <- TRUE
              Other[comp] <- FALSE
              NonClassified[i] <- FALSE
            } else if (3 * df < dt || (df < dt && df > 0.25)) {
              # Unisce il componente a "Forb"
              Forb[comp] <- TRUE
              Other[comp] <- FALSE
              NonClassified[i] <- FALSE
            } else if ((df == 1000 && dt == 700) || dt > k * inputs$PatchDiam1) {
              # Componente isolato, nessuna azione (riprovato al prossimo round)
            } else {
              # Unisce al Trunk
              I_conn <- comp[IC]
              J_conn <- trunk_set[IT]
              Other[comp] <- FALSE
              Trunk[comp] <- TRUE
              NonClassified[i] <- FALSE
              Nei[[I_conn]] <- c(Nei[[I_conn]], J_conn)
              Nei[[J_conn]] <- c(Nei[[J_conn]], I_conn)
            }
          }
        }
      }
      first_round <- FALSE
      if (sum(NonClassified) < npre) {
        again <- TRUE
        npre <- sum(NonClassified)
      } else {
        again <- FALSE
      }
    }
    k <- k + k0
    Cmin <- 3 * Cmin
  }
  Forb[Base] <- FALSE
  cover$neighbor <- Nei
  
  list(cover = cover, Forb = Forb)
}
