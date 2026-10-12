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

# =============================================================================
# tree_sets.R  (versione corretta)
# BUG FIX 1: rimossa ridefinizione di union() che mascherava base::union
# BUG FIX 2: distances_to_line in define_base_forb_with_ground — output è lista
# BUG FIX 3: <<- sostituito con variabili locali in define_trunk
# BUG FIX 4: connected_components ora è in helpers.R
# =============================================================================

#' Determines tree sets, base of trunk, and updates neighbor relations
tree_sets <- function(P, cover, inputs, segment = NULL) {
  np <- nrow(P)
  nb <- length(cover$center)
  
  if (nb == 0) stop("cover_sets produced no cover sets.")
  
  cen <- cover$center
  Ce  <- P[cen, , drop = FALSE]
  nb  <- length(cen)
  
  aux <- list(
    nb     = nb,
    Fal    = rep(FALSE, nb),
    Ind    = seq_len(nb),
    Ce     = Ce,
    Hmin   = min(Ce[, 3]),
    Height = max(Ce[, 3]) - min(Ce[, 3])
  )
  
  if (is.null(segment)) {
    result <- define_base_forb(P, cover, aux, inputs)
  } else {
    inputs$OnlyTree <- TRUE
    result <- define_base_forb(P, cover, aux, inputs, segment)
  }
  Base  <- result$Base
  Forb  <- result$Forb
  cover <- result$cover
  
  if (is.null(segment)) {
    result2 <- define_trunk(cover, aux, Base, Forb, inputs)
  } else {
    result2 <- define_main_branches(cover, segment, aux, inputs)
  }
  Trunk <- result2$Trunk
  cover <- result2$cover
  
  result3 <- make_tree_connected(cover, aux, Forb, Base, Trunk, inputs)
  cover   <- result3$cover
  Forb    <- result3$Forb
  
  list(cover = cover, Base = Base, Forb = Forb)
}

# -----------------------------------------------------------------------------
define_base_forb <- function(P, cover, aux, inputs, segment = NULL) {
  Ce <- aux$Ce
  nb <- aux$nb
  
  if (inputs$OnlyTree && is.null(segment)) {
    BaseHeight <- min(1.5, 0.02 * aux$Height)
    I    <- Ce[, 3] < aux$Hmin + BaseHeight
    Base <- aux$Ind[I]
    
    # FIX: drop=FALSE evita che Ce[Base,1:2] diventi vettore con Base di 1 elemento
    Wb <- if (length(Base) > 1) diff(range(Ce[Base, 1:2, drop=FALSE])) else 0
    Wt <- diff(range(Ce[, 1:2, drop=FALSE]))
    k  <- 1L
    while (k <= 5L && length(Base) > 0 && Wb > 0.3 * Wt) {
      BaseHeight <- max(BaseHeight - 0.05, 0.05)
      I    <- Ce[, 3] < aux$Hmin + BaseHeight
      Base <- aux$Ind[I]
      if (length(Base) > 1) Wb <- diff(range(Ce[Base, 1:2, drop=FALSE]))
      k <- k + 1L
    }
    Forb <- aux$Fal
    
  } else if (inputs$OnlyTree && !is.null(segment)) {
    SoP <- segment$SegmentOfPoint[cover$center]
    cat(sprintf("  [TS2] SoP seg1=%d other=%d\n",
                sum(!is.na(SoP)&SoP==1), sum(!is.na(SoP)&SoP>1)))
    
    BaseHeight <- min(1.5, max(0.05 * aux$Height, 0.5))
    brad <- inputs$BallRad2 %||% inputs$BallRad1 %||% 0.40
    
    # Passo 1: cover di cover2 il cui centro appartiene al segmento 1 (tronco)
    stem1 <- aux$Ind[!is.na(SoP) & SoP == 1]
    cat(sprintf("  [TS2] stem1 (SoP==1): n=%d", length(stem1)))
    if (length(stem1) > 0)
      cat(sprintf(" Z=[%.3f, %.3f]", min(Ce[stem1,3]), max(Ce[stem1,3])))
    cat("\n")
    
    # Passo 2: tra quelli, prendi solo i più bassi (sotto BaseHeight)
    if (length(stem1) > 0) {
      stem_low <- stem1[Ce[stem1,3] < aux$Hmin + BaseHeight]
      if (length(stem_low) == 0) stem_low <- stem1  # tutti se nessuno sotto soglia
    } else {
      # Fallback: i cover più bassi in assoluto
      stem_low <- aux$Ind[order(Ce[aux$Ind,3])[1:min(30, aux$nb)]]
    }
    
    # Passo 3: tra i bassi, stima centro XY del tronco dai punti contenuti
    base_pts <- unlist(cover$ball[stem_low])
    base_pts <- base_pts[!is.na(base_pts) & base_pts >= 1 & base_pts <= nrow(P)]
    if (length(base_pts) > 5) {
      cx <- median(P[base_pts,1]); cy <- median(P[base_pts,2])
      # Raggio del tronco = 75° percentile distanze (esclude outlier)
      dpts <- sqrt((P[base_pts,1]-cx)^2 + (P[base_pts,2]-cy)^2)
      r_trunk_est <- quantile(dpts, 0.75, names=FALSE)
      r_trunk_est <- max(min(r_trunk_est, brad), 0.05)
    } else {
      cx <- mean(Ce[aux$Ind,1]); cy <- mean(Ce[aux$Ind,2])
      r_trunk_est <- brad * 0.5
    }
    
    # Passo 4: aggiunge cover2 bassi vicini al centro tronco non già in stem_low
    all_low <- aux$Ind[Ce[aux$Ind,3] < aux$Hmin + BaseHeight]
    dxy     <- sqrt((Ce[all_low,1]-cx)^2 + (Ce[all_low,2]-cy)^2)
    stem    <- all_low[dxy < 2.0 * r_trunk_est]
    if (length(stem) < 3) stem <- stem_low
    
    # Passo 5: limita la base ai cover più bassi e più centrali (max 5)
    # Una base piccola fa crescere il BFS verticalmente come un filo
    # invece di allargarsi orizzontalmente
    if (length(stem) > 5) {
      # Score = Z bassa + distanza XY piccola dal centro
      dxy_stem <- sqrt((Ce[stem,1]-cx)^2 + (Ce[stem,2]-cy)^2)
      z_stem   <- Ce[stem,3] - aux$Hmin
      # Normalizza su range [0,1] e combina: vogliamo Z bassa e XY piccola
      z_norm  <- if (diff(range(z_stem))   > 1e-6) (z_stem   - min(z_stem))   / diff(range(z_stem))   else rep(0, length(stem))
      xy_norm <- if (diff(range(dxy_stem)) > 1e-6) (dxy_stem - min(dxy_stem)) / diff(range(dxy_stem)) else rep(0, length(stem))
      score   <- z_norm + xy_norm  # minimizza entrambi
      stem    <- stem[order(score)[1:5]]
    }
    
    Base <- stem
    if (length(Base) == 0) Base <- aux$Ind[which.min(Ce[aux$Ind,3])]
    Forb <- aux$Fal
    cat(sprintf("  [TS2] r_trunk_est=%.3f cx=%.2f cy=%.2f\n", r_trunk_est, cx, cy))
    cat(sprintf("  [TS2] stem=%d Base=%d Z=[%.3f, %.3f]\n",
                length(stem), length(Base),
                min(Ce[Base,3]), max(Ce[Base,3])))
    
  } else {
    result <- define_base_forb_with_ground(P, cover, aux, inputs)
    return(result)
  }
  
  list(Base = Base, Forb = Forb, cover = cover)
}

# -----------------------------------------------------------------------------
define_base_forb_with_ground <- function(P, cover, aux, inputs) {
  Ce  <- aux$Ce
  Nei <- cover$neighbor
  Bal <- cover$ball
  
  Min <- apply(Ce, 2, min)
  E   <- min(0.1, 0.015 * aux$Height)
  n_val <- pmax(2L, ceiling((apply(Ce, 2, max)[1:2] - Min[1:2]) / E) + 1L)
  n <- as.integer(n_val)
  
  px <- pmax(1L, pmin(n[1], floor((Ce[, 1] - Min[1]) / E) + 1L))
  py <- pmax(1L, pmin(n[2], floor((Ce[, 2] - Min[2]) / E) + 1L))
  
  LexOrd <- (px - 1L) + (py - 1L) * n[1]
  SortOrd <- order(LexOrd)
  LexOrd  <- LexOrd[SortOrd]
  
  hei <- matrix(0, n[1], n[2])
  den <- hei; baseden <- hei
  Partition <- array(vector("list", prod(n)), dim = n)
  
  p <- 1L
  while (p <= aux$nb) {
    t <- 1L
    while (p + t <= aux$nb && LexOrd[p] == LexOrd[p + t]) t <- t + 1L
    q <- SortOrd[p]
    J <- SortOrd[p:(p + t - 1L)]
    Partition[[px[q], py[q]]] <- list(J)
    K   <- ceiling(10 * (Ce[J, 3] - Min[3] + 0.01) / max(aux$Height - 0.01, 0.01))
    B   <- K <= 2
    K   <- unique(K)
    hei[px[q], py[q]]    <- length(K) / 10
    den[px[q], py[q]]    <- t
    baseden[px[q], py[q]] <- sum(B)
    p <- p + t
  }
  
  den     <- den / max(den, 1)
  baseden <- baseden / max(baseden, 1)
  f       <- den * hei * baseden
  
  for (i in 2:(n[1] - 1)) {
    for (j in 2:(n[2] - 1)) {
      f[i, j] <- mean(f[(i-1):(i+1), (j-1):(j+1)])
    }
  }
  f <- f / max(f, 1)
  
  I      <- f > 0.5
  Trunk0 <- unlist(Partition[I])
  if (length(Trunk0) == 0) Trunk0 <- seq_len(aux$nb)
  
  HBottom <- min(Ce[Trunk0, 3])
  I2 <- Ce[Trunk0, 3] > HBottom + min(0.02 * aux$Height, 0.3) &
    Ce[Trunk0, 3] < HBottom + min(0.08 * aux$Height, 1.5)
  Trunk <- Trunk0[I2]
  if (length(Trunk) == 0) Trunk <- Trunk0
  
  Trunk <- unique(c(Trunk, unlist(Nei[Trunk]), unlist(Nei[unlist(Nei[Trunk])])))
  Trunk <- Trunk[Trunk >= 1 & Trunk <= aux$nb]
  
  # Componenti connesse nel sottografo di Trunk
  {
    nb_loc <- length(Trunk)
    if (nb_loc > 1) {
      idx_map <- integer(aux$nb); idx_map[Trunk] <- seq_len(nb_loc)
      comp_id <- seq_len(nb_loc)
      changed <- TRUE
      while (changed) {
        changed <- FALSE
        for (ii in seq_len(nb_loc)) {
          nb_ii <- Nei[[Trunk[ii]]]
          nb_ii <- nb_ii[!is.na(nb_ii)]
          nb_loc_ii <- idx_map[nb_ii[nb_ii >= 1 & nb_ii <= aux$nb]]
          nb_loc_ii <- nb_loc_ii[nb_loc_ii > 0]
          if (length(nb_loc_ii) == 0) next
          min_c <- min(comp_id[c(ii, nb_loc_ii)])
          if (any(comp_id[c(ii, nb_loc_ii)] != min_c)) {
            comp_id[c(ii, nb_loc_ii)] <- min_c; changed <- TRUE
          }
        }
      }
      Comp <- lapply(unique(comp_id), function(cc) Trunk[comp_id == cc])
      CS   <- sapply(Comp, length)
      Trunk <- Comp[[which.max(CS)]]
    }
  }
  
  I3    <- Ce[Trunk, 3] < HBottom + min(0.1 * aux$Height, 2)
  Trunk <- Trunk[I3]
  Trunk <- unique(c(Trunk, unlist(Nei[Trunk])))
  Trunk <- Trunk[Trunk >= 1 & Trunk <= aux$nb]
  
  # FIX: distances_to_line ritorna lista — usa $d
  Points  <- Ce[Trunk, , drop = FALSE]
  cstart  <- colMeans(Points)
  dl      <- distances_to_line(Points, c(0, 0, 1), cstart)
  cradius <- mean(dl$d)
  
  BaseHeight <- min(1.5, 0.02 * aux$Height)
  Bot  <- min(Ce[Trunk, 3])
  J    <- Ce[Trunk, 3] < Bot + BaseHeight
  Base <- Trunk[J]
  if (length(Base) == 0) Base <- Trunk[which.min(Ce[Trunk, 3])]
  
  Trunk2 <- unique(c(Trunk, unlist(Nei[Trunk])))
  Trunk2 <- Trunk2[Trunk2 >= 1 & Trunk2 <= aux$nb]
  Forb   <- aux$Fal
  Ground <- setdiff(unlist(Nei[Base]), Trunk2)
  Ground <- setdiff(base::union(Ground, unlist(Nei[Ground])), Trunk2)
  Forb[Ground] <- TRUE
  Forb[Base]   <- FALSE
  
  Forb[Ce[, 3] < Bot + 0.03 * aux$Height] <- TRUE
  Forb[Trunk2] <- FALSE
  
  cover$ball <- Bal
  list(Base = Base, Forb = Forb, cover = cover)
}

# -----------------------------------------------------------------------------
define_trunk <- function(cover, aux, Base, Forb, inputs) {
  Nei <- cover$neighbor
  Ce  <- aux$Ce
  
  Trunk <- aux$Fal
  Trunk[Base] <- TRUE
  
  L    <- if (!is.null(inputs$PatchDiam1)) inputs$PatchDiam1 else 0.1
  H0   <- min(Ce[Base, 3], na.rm = TRUE)
  H    <- max(Ce[Base, 3], na.rm = TRUE)
  Exp  <- Base
  FirstMod <- TRUE
  
  # Variabili locali per distanza/coseno (FIX: era <<-)
  Dist_mat <- matrix(3, 1, 1)
  Cos_mat  <- matrix(0, 1, 1)
  
  repeat {
    Exp0 <- Exp
    Exp  <- unique(c(Exp, unlist(Nei[Exp])))
    Exp  <- Exp[!Trunk[Exp] & !Forb[Exp]]
    
    if (length(Exp) > 0) {
      Trunk[Exp] <- TRUE
      H <- max(Ce[Exp, 3], na.rm = TRUE)
    }
    
    stop_exp    <- (length(Exp) == 0)
    height_stuck <- !is.na(H) && !is.na(H0) && (H < H0 + L / 2)
    below_canopy <- !is.na(H) && (H < aux$Hmin + aux$Height - 1)
    
    if ((stop_exp || height_stuck) && below_canopy) {
      if (FirstMod) {
        FirstMod <- FALSE
      }
      Region     <- if (length(Exp) > 0) Exp else Exp0
      if (is.na(H)) break
      
      SetsAbove <- which(Ce[, 3] > H & Ce[, 3] < H + 2 & !Trunk & !Forb)
      
      if (length(SetsAbove) > 0 && length(Region) > 0) {
        n_inner <- length(Region)
        m_inner <- length(SetsAbove)
        Dist_mat <- matrix(3, n_inner, m_inner)
        Cos_mat  <- matrix(0, n_inner, m_inner)
        
        for (i in seq_len(n_inner)) {
          V_diff  <- sweep(Ce[SetsAbove, , drop = FALSE], 2, Ce[Region[i], ])
          LenSq   <- rowSums(V_diff^2)
          LenSq[LenSq < 1e-10] <- 1e-10
          Dist_mat[i, ] <- LenSq
          Cos_mat[i, ]  <- V_diff[, 3] / sqrt(LenSq)
        }
        
        I_sel <- Cos_mat > 0.7
        t_count <- 0L
        while (!any(I_sel) && t_count < 10L) {
          t_count <- t_count + 1L
          I_sel <- Cos_mat > (0.7 - t_count * 0.05)
        }
        
        if (any(I_sel)) {
          Dist_mat[!I_sel] <- 3
          J_idx <- which.min(apply(Dist_mat, 1, min))
          I_idx <- which.min(Dist_mat[J_idx, ])
          
          t_node <- Region[J_idx]
          s_node <- SetsAbove[I_idx]
          
          if (!is.na(t_node) && !is.na(s_node)) {
            Nei[[t_node]] <- unique(c(Nei[[t_node]], s_node))
            Nei[[s_node]] <- unique(c(Nei[[s_node]], t_node))
            Exp    <- s_node
            Trunk[Exp] <- TRUE
            H <- max(Ce[Exp, 3], na.rm = TRUE)
          } else break
        } else break
      } else break
      
    } else if (stop_exp) {
      break
    }
    
    if (!is.na(H)) H0 <- H
  }
  
  cover$neighbor <- Nei
  list(Trunk = Trunk, cover = cover)
}

# -----------------------------------------------------------------------------
define_main_branches <- function(cover, segment, aux, inputs) {
  Bal <- cover$ball
  Nei <- cover$neighbor
  nb  <- length(Bal)
  
  SoP     <- segment$SegmentOfPoint
  max_seg <- max(SoP, na.rm = TRUE)
  if (max_seg < 1) {
    return(list(Trunk = aux$Fal, cover = cover))
  }
  
  MBI <- rep(FALSE, max_seg)
  MBI[1] <- TRUE
  for (idx in c(segment$branch1indexes, segment$branch2indexes, segment$branch3indexes)) {
    if (!is.na(idx) && idx >= 1 && idx <= max_seg) MBI[idx] <- TRUE
  }
  
  MainBranches <- rep(FALSE, nb)
  for (i in seq_len(nb)) {
    pts <- Bal[[i]]
    if (is.null(pts) || length(pts) == 0) next
    bi  <- unique(SoP[pts])
    bi  <- bi[!is.na(bi) & bi >= 1 & bi <= max_seg]
    if (length(bi) > 0 && any(MBI[bi])) MainBranches[i] <- TRUE
  }
  
  Trunk <- aux$Fal
  Trunk[MainBranches] <- TRUE
  cover$neighbor <- Nei
  list(Trunk = Trunk, cover = cover)
}

# -----------------------------------------------------------------------------
make_tree_connected <- function(cover, aux, Forb, Base, Trunk, inputs) {
  Nei <- cover$neighbor
  
  Trunk[Forb] <- FALSE
  Exp <- which(Trunk)
  
  # FIX: il loop originale usava vettori logici con indicizzazione mista
  changed <- TRUE
  while (changed) {
    new_trunk <- unique(unlist(Nei[Exp]))
    new_trunk <- new_trunk[new_trunk >= 1 & new_trunk <= aux$nb]
    new_trunk <- new_trunk[!Trunk[new_trunk] & !Forb[new_trunk]]
    if (length(new_trunk) == 0) { changed <- FALSE; break }
    Trunk[new_trunk] <- TRUE
    Exp <- new_trunk
  }
  
  Other <- !Forb & !Trunk
  Other[Base] <- FALSE
  
  if (any(Other)) {
    # FIX: marca come Forb i cover non connessi all'albero
    # Componenti connesse nel sottografo Other: marca non-connesse come Forb
    {
      other_nodes <- which(Other)
      nb_o <- length(other_nodes)
      if (nb_o > 0) {
        idx_o <- integer(aux$nb); idx_o[other_nodes] <- seq_len(nb_o)
        comp_o <- seq_len(nb_o)
        changed_o <- TRUE
        while (changed_o) {
          changed_o <- FALSE
          for (ii in seq_len(nb_o)) {
            nb_ii <- Nei[[other_nodes[ii]]]
            nb_ii <- nb_ii[!is.na(nb_ii)]
            nb_loc_ii <- idx_o[nb_ii[nb_ii >= 1 & nb_ii <= aux$nb]]
            nb_loc_ii <- nb_loc_ii[nb_loc_ii > 0]
            if (length(nb_loc_ii) == 0) next
            min_c <- min(comp_o[c(ii, nb_loc_ii)])
            if (any(comp_o[c(ii, nb_loc_ii)] != min_c)) {
              comp_o[c(ii, nb_loc_ii)] <- min_c; changed_o <- TRUE
            }
          }
        }
        Comps <- lapply(unique(comp_o), function(cc) other_nodes[comp_o == cc])
        tree_comp <- which(sapply(Comps, function(comp) any(comp %in% which(Trunk))))
        if (length(tree_comp) > 0) {
          tree_set <- unlist(Comps[tree_comp])
          non_tree <- setdiff(other_nodes, tree_set)
          if (length(non_tree) > 0) Forb[non_tree] <- TRUE
        }
      }
    }
  }
  
  Forb[Base] <- FALSE
  cover$neighbor <- Nei
  list(cover = cover, Forb = Forb)
}
