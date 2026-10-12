# =============================================================================
# CORREZIONI SESSIONE 2026 - Day 20/21
# =============================================================================
# Questo file raccoglie tutte le correzioni scoperte e verificate in questa
# sessione, che altrimenti esisterebbero solo in memoria e andrebbero perse
# al riavvio di R. NON modifica i file originali (piu' sicuro e tracciabile):
# va caricato con source() PER ULTIMO, dopo tutti gli altri file del progetto
# (incluso pipeline_integrata_finale.R), cosi' le sue definizioni vincono.
#
# ESCLUSIONI DAL CARICAMENTO (lista da riusare in ogni sessione, non codice
# da caricare qui): tools.R, extract_trees.R, estrai_alberi_qsm.R, le tre
# versioni scartate (cover_sets_MYFIX_backup.R, cover_sets_original.R,
# treeqsm_patched.R, identify_stem_chain(1).R), gli script di test/debug,
# i quattro file pipeline_*.R di alto livello (SOSTITUITI da questo file
# insieme a pipeline_integrata_finale.R).
# =============================================================================

# --- FIX 1: gruppo coerente di funzioni geometriche da helpers.R -----------
# least_squares_fitting.R (file satellite, mai risolto il conflitto con
# helpers.R per ordine alfabetico) manca delle derivate DR1/DR2 in
# form_rotation_matrices(), necessarie al Jacobiano di func_grad_cylinder().
# Verificato sperimentalmente con un cilindro sintetico di raggio noto
# (0.15m): la versione di least_squares_fitting.R stima 0.53 (conv=FALSE),
# questa versione stima 0.15 esatto (conv=TRUE).

rotation_matrix <- function(A, angle) {
  A <- A / norm_vec(A)
  c <- cos(angle); s <- sin(angle)
  R <- matrix(0, 3, 3)
  R[1, ] <- c(A[1]^2 + (1 - A[1]^2) * c, A[1] * A[2] * (1 - c) - A[3] * s, A[1] * A[3] * (1 - c) + A[2] * s)
  R[2, ] <- c(A[1] * A[2] * (1 - c) + A[3] * s, A[2]^2 + (1 - A[2]^2) * c, A[2] * A[3] * (1 - c) - A[1] * s)
  R[3, ] <- c(A[1] * A[3] * (1 - c) - A[2] * s, A[2] * A[3] * (1 - c) + A[1] * s, A[3]^2 + (1 - A[3]^2) * c)
  R
}

rotate_to_z_axis <- function(Vec) {
  D <- cross_product(Vec, c(0, 0, 1))
  if (norm_vec(D) > 0) {
    a <- acos(Vec[3]); R <- rotation_matrix(D, a)
  } else { R <- diag(3); a <- 0; D <- c(1, 0, 0) }
  list(R = R, D = D, a = a)
}

form_rotation_matrices <- function(angles) {
  alpha <- angles[1]; beta <- angles[2]
  ca <- cos(alpha); sa <- sin(alpha); cb <- cos(beta); sb <- sin(beta)
  R <- matrix(0, 3, 3)
  R[1, ] <- c(cb, 0, sb); R[2, ] <- c(sa * sb, ca, -sa * cb); R[3, ] <- c(-ca * sb, sa, ca * cb)
  DR1 <- matrix(0, 3, 3); DR1[2, ] <- c(ca * sb, -sa, -ca * cb); DR1[3, ] <- c(sa * sb, ca, sa * cb)
  DR2 <- matrix(0, 3, 3); DR2[1, ] <- c(-sb, 0, cb); DR2[2, ] <- c(sa * cb, 0, sa * sb); DR2[3, ] <- c(-ca * cb, 0, -ca * sb)
  list(R = R, DR1 = DR1, DR2 = DR2)
}

func_grad_cylinder <- function(par, P, weight = NULL) {
  x0 <- par[1]; y0 <- par[2]; alpha <- par[3]; beta <- par[4]; r <- par[5]
  rot <- form_rotation_matrices(c(alpha, beta)); R <- rot$R; DR1 <- rot$DR1; DR2 <- rot$DR2
  Pt  <- (P - matrix(c(x0, y0, 0), nrow = nrow(P), ncol = 3, byrow = TRUE)) %*% t(R)
  xt  <- Pt[, 1]; yt <- Pt[, 2]; rt <- sqrt(xt^2 + yt^2); rt_safe <- pmax(rt, 1e-10)
  dist <- rt - r
  if (!is.null(weight)) dist <- weight * dist
  n <- nrow(P); J <- matrix(0, n, 5); N <- cbind(xt / rt_safe, yt / rt_safe); A <- diag(2)
  A1 <- R %*% c(-1, 0, 0); A[1,1] <- A1[1]; A[2,2] <- A1[2]; J[, 1] <- rowSums(N %*% A)
  A2 <- R %*% c(0, -1, 0); A[1,1] <- A2[1]; A[2,2] <- A2[2]; J[, 2] <- rowSums(N %*% A)
  if (!is.null(DR1)) { dPt1 <- (P - matrix(c(x0, y0, 0), nrow = n, ncol = 3, byrow = TRUE)) %*% t(DR1); J[, 3] <- rowSums(N * dPt1[, 1:2, drop = FALSE]) }
  if (!is.null(DR2)) { dPt2 <- (P - matrix(c(x0, y0, 0), nrow = n, ncol = 3, byrow = TRUE)) %*% t(DR2); J[, 4] <- rowSums(N * dPt2[, 1:2, drop = FALSE]) }
  J[, 5] <- -1
  if (!is.null(weight)) J <- J * weight
  list(d = dist, J = J)
}

least_squares_cylinder <- function(P, cyl0, weight = NULL, Q = NULL) {
  iter <- 0L; conv <- FALSE; rel <- TRUE; NoWeights <- is.null(weight)
  rot <- rotate_to_z_axis(cyl0$axis); Rot0 <- rot$R
  Pt <- (P - matrix(cyl0$start, nrow = nrow(P), ncol = 3, byrow = TRUE)) %*% t(Rot0)
  par <- c(1e-6, 1e-6, 1e-6, 1e-6, cyl0$radius); SS_prev <- Inf; dist <- rep(0, nrow(P))
  while (iter < 50L && !conv) {
    fg <- if (NoWeights) func_grad_cylinder(par, Pt) else func_grad_cylinder(par, Pt, weight)
    d0 <- fg$d; J <- fg$J
    if (is.null(J) || any(!is.finite(J))) { rel <- FALSE; break }
    A <- t(J) %*% J; b <- t(J) %*% d0; SS0 <- sqrt(sum(d0^2))
    lam <- 1e-4 * max(diag(A), na.rm = TRUE); A_lm <- A + lam * diag(5)
    p_step <- tryCatch(suppressWarnings(-solve(A_lm, b)), error = function(e) NULL)
    if (is.null(p_step) || any(!is.finite(p_step))) { rel <- FALSE; break }
    par_new <- par + p_step
    dist_new <- if (NoWeights) func_grad_cylinder(par_new, Pt)$d else func_grad_cylinder(par_new, Pt, weight)$d
    SS1 <- sqrt(sum(dist_new^2))
    if (SS1 <= SS0) { par <- par_new; dist <- dist_new } else { dist <- d0 }
    if (rcond(A_lm) < 1e-12) rel <- FALSE
    if (abs(SS_prev - SS1) < 1e-5) conv <- TRUE
    SS_prev <- SS1; iter <- iter + 1L
  }
  rot2 <- form_rotation_matrices(par[3:4])
  Axis <- as.numeric(t(Rot0) %*% t(rot2$R) %*% c(0, 0, 1))
  Point <- as.numeric(t(Rot0) %*% c(par[1], par[2], 0) + cyl0$start)
  P_use <- if (!is.null(Q) && nrow(Q) > 5) Q else P
  H <- as.numeric(P_use %*% Axis); hmin <- min(H)
  Point <- Point - (sum(Axis * Point) - hmin) * Axis
  if (!is.null(weight)) { I_w <- weight == max(weight); mad_val <- mean(abs(dist[I_w])) } else { mad_val <- mean(abs(dist)) }
  cyl <- list(radius = as.numeric(par[5]), length = as.numeric(abs(max(H) - hmin)), start = as.numeric(Point),
              axis = as.numeric(Axis), mad = if (is.finite(mad_val)) mad_val else 0, conv = conv, rel = rel, SurfCov = 0)
  if (is.finite(cyl$radius) && cyl$radius > 1e-4 && is.finite(cyl$length) && cyl$length > 1e-4 && !any(is.na(Axis)) && !any(is.na(Point))) {
    res_l <- max(3L, ceiling(cyl$length / 0.03)); res_s <- min(36L, max(8L, ceiling(2 * pi * cyl$radius / 0.03)))
    cyl$SurfCov <- tryCatch(surface_coverage(P_use, Axis, Point, res_l, res_s, 0.8 * cyl$radius), error = function(e) 0)
  }
  cyl
}

# --- FIX 2: adjustments() - divisione per zero quando sum(SC)==0 -----------
# Riga ~160 dell'originale: quando nessun cilindro di un segmento corto
# raggiunge SurfCov>=0.7, "a <- sum(R*SC/sum(SC))" e' 0/0=NaN se tutti i
# pesi sono zero, rompendo seq(a,b,...) a valle. Confermato riproducendo
# l'errore con un caso sintetico (3 cilindri, tutti SurfCov=0) e verificando
# che il fix lo risolva.

adjustments <- function(cyl, parcyl, inputs, Regs) {
  nc <- length(cyl$radius)
  if (nc == 0) return(cyl)
  Mod <- rep(FALSE, nc)
  SC  <- cyl$SurfCov
  MinCylRad <- inputs$MinCylRad %||% 0.005
  has_parent <- !is.null(parcyl$radius) && length(parcyl$radius) > 0 && is.finite(parcyl$radius[1])
  if (has_parent) {
    MaxR <- 0.95 * parcyl$radius[1]; MaxR <- max(MaxR, MinCylRad)
  } else {
    a <- min(3, nc); MaxR <- 1.25 * max(cyl$radius[seq_len(a)])
  }
  MinR <- if (any(SC > 0.7)) min(cyl$radius[SC > 0.7]) else NA_real_
  if (!is.na(MinR) && min(cyl$radius) < MinR / 2) {
    MinR <- if (any(SC > 0.4)) min(cyl$radius[SC > 0.4]) else NA_real_
  } else if (is.na(MinR)) {
    MinR <- if (any(SC > 0.4)) min(cyl$radius[SC > 0.4]) else MinCylRad
  }
  if (is.na(MinR)) MinR <- MinCylRad
  I <- cyl$radius < MinR; cyl$radius[I] <- MinR; Mod[I] <- TRUE
  ParentCor <- isTRUE(inputs$ParentCor %||% TRUE)
  if (ParentCor || nc <= 3) {
    I <- (cyl$radius > MaxR & SC < 0.7) | (cyl$radius > 1.2 * MaxR)
    cyl$radius[I] <- MaxR; Mod[I] <- TRUE
    if (nc <= 3) {
      I <- (cyl$radius > 0.75 * MaxR & SC < 0.7)
      if (any(I)) { r <- pmax(SC[I] / 0.7 * cyl$radius[I], MinR); cyl$radius[I] <- r; Mod[I] <- TRUE }
    }
  }
  if (nc > 2) {
    for (i in 2:(nc - 1)) {
      if (SC[i] < 0.7 && SC[i - 1] >= 0.7 && SC[i + 1] >= 0.7) {
        cyl$radius[i] <- 0.5 * (cyl$radius[i - 1] + cyl$radius[i + 1]); Mod[i] <- TRUE
      }
    }
  }
  TaperCor <- isTRUE(inputs$TaperCor %||% TRUE)
  if (TaperCor) {
    if (max(cyl$radius) < 0.001) {
      if (nc > 2) {
        r <- sort(cyl$radius); r <- r[2:(length(r) - 1)]; a <- 2 * mean(r)
        if (a > max(r)) a <- min(0.01, max(r))
        b <- min(0.5 * min(cyl$radius), 0.001)
        cyl$radius <- seq(a, b, length.out = nc)
      } else if (nc > 1) {
        r <- max(cyl$radius); cyl$radius <- c(r, 0.5 * r)
      }
      Mod <- rep(TRUE, nc)
    } else if (nc > 4) {
      branchlen <- sum(cyl$length[seq_len(nc)])
      L <- cyl$length / 2 + c(0, cumsum(cyl$length[seq_len(nc - 1)]))
      Taper1 <- c(L, branchlen); Taper2 <- c(1.05 * cyl$radius, MinR); sc <- c(SC, 1)
      A_mat <- matrix(c(sum(sc * Taper1^4), sum(sc * Taper1^2), sum(sc * Taper1^2), sum(sc)), 2, 2, byrow = TRUE)
      y_vec <- c(sum(sc * Taper2 * Taper1^2), sum(sc * Taper2))
      x_sol <- tryCatch(solve(A_mat, y_vec), error = function(e) c(-0.0001, MinR))
      x_sol[1] <- min(x_sol[1], -0.0001)
      Ru <- x_sol[1] * L^2 + x_sol[2]; Ru[Ru < MinR] <- MinR
      if (max(Ru) > MaxR) { a <- max(Ru); Ru <- MaxR / a * Ru }
      Rl <- 0.75 * Ru; Rl[Rl < MinR] <- MinR
      I <- cyl$radius > Ru & SC < 0.7; cyl$radius[I] <- Ru[I] + (cyl$radius[I] - Ru[I]) * SC[I] / 0.7; Mod[I] <- TRUE
      I <- cyl$radius > 1.333 * Ru & SC >= 0.7; cyl$radius[I] <- Ru[I] + (cyl$radius[I] - Ru[I]) * SC[I]; Mod[I] <- TRUE
      I <- (cyl$radius < Rl & SC < 0.7) | (cyl$radius < 0.5 * Rl); cyl$radius[I] <- Rl[I]; Mod[I] <- TRUE
    } else {
      R <- cyl$radius
      if (sum(SC >= 0.7) > 1) {
        a <- max(R[SC >= 0.7]); b <- min(R[SC >= 0.7])
      } else if (sum(SC >= 0.7) == 1) {
        a <- max(R[SC >= 0.7]); b <- min(R)
      } else {
        a <- if (sum(SC) > 0) sum(R * SC / sum(SC)) else mean(R)
        b <- min(R)
      }
      Ru <- seq(a, b, length.out = nc)
      I <- SC < 0.7 & !Mod; cyl$radius[I] <- Ru[I] + (R[I] - Ru[I]) * SC[I] / 0.7; Mod[I] <- TRUE
    }
  }
  nr <- length(Regs); Reg <- NULL
  for (i in seq_len(nc)) {
    if (Mod[i]) {
      if (nr == nc) { Reg <- Regs[[i]] } else if (i > 1) { Reg <- Regs[[i - 1]] }
      if (!is.null(Reg) && abs(cyl$radius[i] - cyl$radius0[i]) > 0.005 && (nr == nc || (nr < nc && i > 1))) {
        P <- sweep(Reg, 2, cyl$start[i, ]); ow <- orthonormal_vectors(cyl$axis[i, ]); U <- ow$V; V <- ow$W
        P_proj <- cbind(as.vector(P %*% U), as.vector(P %*% V))
        cir <- tryCatch(least_squares_circle_centre(P_proj, c(0, 0), cyl$radius[i]), error = function(e) NULL)
        if (!is.null(cir) && isTRUE(cir$conv) && isTRUE(cir$rel)) {
          cyl$start[i, ] <- cyl$start[i, ] + cir$point[1] * U + cir$point[2] * V
          cyl$mad[i] <- cir$mad
          dl <- distances_to_line(Reg, cyl$axis[i, ], cyl$start[i, ]); h <- dl$h
          if (min(h) < -0.001) {
            cyl$length[i] <- max(h) - min(h)
            cyl$start[i, ] <- cyl$start[i, ] + min(h) * cyl$axis[i, ]
            dl <- distances_to_line(Reg, cyl$axis[i, ], cyl$start[i, ]); h <- dl$h
          }
          a <- max(0.02, 0.2 * cyl$radius[i]); nl <- ceiling(cyl$length[i] / a); nl <- max(nl, 4)
          ns <- ceiling(2 * pi * cyl$radius[i] / a); ns <- max(ns, 10); ns <- min(ns, 36)
          cyl$SurfCov[i] <- surface_coverage2(cyl$axis[i, ], cyl$length[i], dl$V, h, nl, ns)
        }
      }
    }
  }
  if (nc > 1) {
    for (j in 2:nc) {
      U <- cyl$start[j, ] - cyl$start[j - 1, ] - cyl$length[j - 1] * cyl$axis[j - 1, ]
      if (sqrt(sum(U^2)) > 0.0001) {
        N <- cyl$axis[j, ]
        if (sqrt(sum(N^2)) > 0) {
          ow <- orthonormal_vectors(N); Vv <- ow$V; Ww <- ow$W; M <- cbind(N, Vv, Ww)
          x <- tryCatch(solve(M, U), error = function(e) c(0, 0, 0))
          cyl$start[j, ] <- cyl$start[j, ] - x[1] * N
          if (x[1] > 0) { cyl$length[j] <- cyl$length[j] + x[1] }
          else if (cyl$length[j] + x[1] > 0) { cyl$length[j] <- cyl$length[j] + x[1] }
        }
      }
    }
  }
  if (has_parent) {
    dl <- distances_to_line(matrix(cyl$start[1, ], 1, 3), parcyl$axis[1, ], parcyl$start[1, ])
    d <- dl$d[1] - parcyl$radius[1]
    if (is.finite(d) && d > 0.001) {
      Vn <- dl$V[1, ]; Vn <- parcyl$radius[1] * Vn / sqrt(sum(Vn^2)); h <- dl$h[1]; B <- dl$B[1, ]
      taper <- cyl$start[1, ]; E <- taper + cyl$length[1] * cyl$axis[1, ]
      if (h >= 0 && h <= parcyl$length[1]) { cyl$start[1, ] <- parcyl$start[1, ] + B + Vn }
      else if (h < 0) { cyl$start[1, ] <- parcyl$start[1, ] + Vn }
      else { cyl$start[1, ] <- parcyl$start[1, ] + parcyl$length[1] * parcyl$axis[1, ] + Vn }
      new_axis <- E - cyl$start[1, ]; new_len <- sqrt(sum(new_axis^2))
      if (new_len > 1e-9) { cyl$axis[1, ] <- new_axis / new_len; cyl$length[1] <- new_len }
    }
  }
  cyl
}

# --- FIX 3: 4 funzioni duplicate tra tree_data_distributions.R (protetta)
# e update_tree_data.R (fragile, vince per ordine alfabetico: "u" dopo "t").
# La versione fragile manca na.rm=TRUE in max(100*Par) e usa "for(i in 1:n)"
# invece di "for(i in seq_len(n))" - se n=0 o Par contiene NA, propaga NA
# fino a matrix(0,8,NA), causando i warning sistematici osservati.

cylinder_distribution <- function(treedata, cyl, dist) {
  if (dist == "Dia") {
    Par <- cyl$radius; n <- ceiling(max(200 * cyl$radius, na.rm = TRUE)); a <- 0.005
  } else if (dist == "Zen") {
    Par <- 180 / pi * acos(pmax(-1, pmin(1, cyl$axis[, 3]))); n <- 18; a <- 10
  } else if (dist == "Azi") {
    Par <- 180 / pi * atan2(cyl$axis[, 2], cyl$axis[, 1]) + 180; n <- 36; a <- 10
  } else { stop("cylinder_distribution: dist non riconosciuto") }
  if (!is.finite(n) || n < 1) n <- 1
  CylDist <- matrix(0, 3, n)
  for (i in seq_len(n)) {
    K <- !is.na(Par) & Par >= (i - 1) * a & Par < i * a
    CylDist[1, i] <- 1000 * pi * sum(cyl$radius[K]^2 * cyl$length[K])
    CylDist[2, i] <- 2 * pi * sum(cyl$radius[K] * cyl$length[K])
    CylDist[3, i] <- sum(cyl$length[K])
  }
  treedata[[paste0("VolCyl", dist)]] <- CylDist[1, ]; treedata[[paste0("AreCyl", dist)]] <- CylDist[2, ]; treedata[[paste0("LenCyl", dist)]] <- CylDist[3, ]
  treedata
}

cylinder_height_distribution <- function(treedata, cylinder, ind) {
  Rad <- cylinder$radius; Len <- cylinder$length; Axe <- cylinder$axis
  MaxHei <- ceiling(treedata$TreeHeight)
  if (!is.finite(MaxHei) || MaxHei < 1) MaxHei <- 1
  treedata$VolCylHei <- numeric(MaxHei); treedata$AreCylHei <- numeric(MaxHei); treedata$LenCylHei <- numeric(MaxHei)
  End <- cylinder$start + cbind(Len * Axe[, 1], Len * Axe[, 2], Len * Axe[, 3])
  bot <- min(cylinder$start[, 3])
  B <- cylinder$start[, 3] - bot; Top <- End[, 3] - bot
  for (j in seq_len(MaxHei)) {
    I1 <- B >= (j - 2) & B < (j - 1); J1 <- B >= (j - 1) & B < j; K1 <- B >= j & B < (j + 1)
    I2 <- Top >= (j - 2) & Top < (j - 1); J2 <- Top >= (j - 1) & Top < j; K2 <- Top >= j & Top < (j + 1)
    C1 <- ind[J1 & J2]; C2 <- ind[J1 & K2]; C3 <- ind[J1 & I2]; C4 <- ind[I1 & J2]; C5 <- ind[K1 & J2]
    v1 <- 1000 * pi * sum(Rad[C1]^2 * Len[C1]); a1 <- 2 * pi * sum(Rad[C1] * Len[C1]); l1 <- sum(Len[C1])
    r2 <- (j - B[C2]) / (Top[C2] - B[C2]); v2 <- 1000 * pi * sum(Rad[C2]^2 * Len[C2] * r2); a2 <- 2 * pi * sum(Rad[C2] * Len[C2] * r2); l2 <- sum(Len[C2] * r2)
    r3 <- (B[C3] - j + 1) / (B[C3] - Top[C3]); v3 <- 1000 * pi * sum(Rad[C3]^2 * Len[C3] * r3); a3 <- 2 * pi * sum(Rad[C3] * Len[C3] * r3); l3 <- sum(Len[C3] * r3)
    r4 <- (Top[C4] - j + 1) / (Top[C4] - B[C4]); v4 <- 1000 * pi * sum(Rad[C4]^2 * Len[C4] * r4); a4 <- 2 * pi * sum(Rad[C4] * Len[C4] * r4); l4 <- sum(Len[C4] * r4)
    r5 <- (j - Top[C5]) / (B[C5] - Top[C5]); v5 <- 1000 * pi * sum(Rad[C5]^2 * Len[C5] * r5); a5 <- 2 * pi * sum(Rad[C5] * Len[C5] * r5); l5 <- sum(Len[C5] * r5)
    treedata$VolCylHei[j] <- v1+v2+v3+v4+v5; treedata$AreCylHei[j] <- a1+a2+a3+a4+a5; treedata$LenCylHei[j] <- l1+l2+l3+l4+l5
  }
  treedata
}

branch_distribution <- function(treedata, branch, dist) {
  n_b <- length(branch$order)
  empty_out <- function(td) {
    for (nm in c("VolBranch", "VolBranch1", "AreBranch", "AreBranch1", "LenBranch", "LenBranch1", "NumBranch", "NumBranch1")) {
      td[[paste0(nm, dist)]] <- numeric(0)
    }
    td
  }
  if (n_b < 2) return(empty_out(treedata))
  BOrd <- branch$order[2:n_b]; BVol <- branch$volume[2:n_b]; BAre <- branch$area[2:n_b]; BLen <- branch$length[2:n_b]
  if (dist == "Dia") { Par <- branch$diameter[2:n_b]; n <- ceiling(max(100 * Par, na.rm = TRUE)); a <- 0.005 }
  else if (dist == "Hei") { Par <- branch$height[2:n_b]; n <- ceiling(treedata$TreeHeight); a <- 1 }
  else if (dist == "Ang") { Par <- branch$angle[2:n_b]; n <- 18; a <- 10 }
  else if (dist == "Zen") { Par <- branch$zenith[2:n_b]; n <- 18; a <- 10 }
  else if (dist == "Azi") { Par <- branch$azimuth[2:n_b] + 180; n <- 36; a <- 10 }
  else { stop("branch_distribution: dist non riconosciuto") }
  if (!is.finite(n) || is.na(n)) n <- 0
  n <- max(0, n)
  BranchDist <- matrix(0, 8, max(n, 1))
  if (n > 0) {
    for (i in seq_len(n)) {
      I <- !is.na(Par) & Par >= (i - 1) * a & Par < i * a
      I1 <- I & !is.na(BOrd) & BOrd == 1
      BranchDist[1, i] <- sum(BVol[I]); BranchDist[2, i] <- sum(BVol[I1])
      BranchDist[3, i] <- sum(BAre[I]); BranchDist[4, i] <- sum(BAre[I1])
      BranchDist[5, i] <- sum(BLen[I]); BranchDist[6, i] <- sum(BLen[I1])
      BranchDist[7, i] <- sum(I); BranchDist[8, i] <- sum(I1)
    }
    BranchDist <- BranchDist[, seq_len(n), drop = FALSE]
  } else { BranchDist <- matrix(0, 8, 0) }
  treedata[[paste0("VolBranch", dist)]]  <- BranchDist[1, ]; treedata[[paste0("VolBranch1", dist)]] <- BranchDist[2, ]
  treedata[[paste0("AreBranch", dist)]]  <- BranchDist[3, ]; treedata[[paste0("AreBranch1", dist)]] <- BranchDist[4, ]
  treedata[[paste0("LenBranch", dist)]]  <- BranchDist[5, ]; treedata[[paste0("LenBranch1", dist)]] <- BranchDist[6, ]
  treedata[[paste0("NumBranch", dist)]]  <- BranchDist[7, ]; treedata[[paste0("NumBranch1", dist)]] <- BranchDist[8, ]
  treedata
}

branch_order_distribution <- function(treedata, branch) {
  BO <- suppressWarnings(max(branch$order, na.rm = TRUE))
  if (!is.finite(BO)) BO <- 0
  n_rows <- max(1, BO)
  BranchOrdDist <- matrix(0, n_rows, 4)
  for (i in seq_len(n_rows)) {
    I <- !is.na(branch$order) & branch$order == i
    BranchOrdDist[i, 1] <- sum(branch$volume[I]); BranchOrdDist[i, 2] <- sum(branch$area[I])
    BranchOrdDist[i, 3] <- sum(branch$length[I]); BranchOrdDist[i, 4] <- sum(I)
  }
  treedata$VolBranchOrd <- BranchOrdDist[, 1]; treedata$AreBranchOrd <- BranchOrdDist[, 2]
  treedata$LenBranchOrd <- BranchOrdDist[, 3]; treedata$NumBranchOrd <- BranchOrdDist[, 4]
  treedata
}

# --- FIX 4: funzioni di supporto costruite in questa sessione --------------

appiattisci_inputs <- function(inputs) {
  inputs$PatchDiam1 <- inputs$PatchDiam1[1]
  inputs$PatchDiam2Min <- inputs$PatchDiam2Min[1]
  inputs$PatchDiam2Max <- inputs$PatchDiam2Max[1]
  inputs$BallRad1 <- inputs$BallRad1[1]
  inputs$BallRad2 <- inputs$BallRad2[1]
  inputs
}

modello_non_degenere <- function(QSM_list, min_height_frac = 0.5, min_cilindri = 3) {
  if (length(QSM_list) == 0) return(NULL)
  non_null <- sapply(QSM_list, function(q) !is.null(q) && !is.null(q$treedata))
  if (!any(non_null)) return(NULL)
  heights <- sapply(QSM_list[non_null], function(q) {
    h <- q$treedata$TreeHeight
    if (is.null(h) || length(h) == 0 || !is.finite(h)) NA_real_ else h
  })
  med_h <- median(heights, na.rm = TRUE)
  ok_height <- !is.na(heights) & (is.na(med_h) | heights >= min_height_frac * med_h)
  n_cyl <- sapply(QSM_list[non_null], function(q) length(q$cylinder$radius))
  ok_cyl <- n_cyl >= min_cilindri
  ok_tutti <- ok_height & ok_cyl
  candidates_idx <- which(non_null)[ok_tutti]
  if (length(candidates_idx) == 0) candidates_idx <- which(non_null)
  QSM_list[[candidates_idx[1]]]
}

# --- Verifiche dopo il riavvio e il bootstrap ---
# mk <- function(d) list(treedata=list(DBHcyl=d/100, TreeHeight=10), cylinder=list(radius=rep(0.1,10)))
# dbh_robusto(lapply(c(1,1,1,1,21.4,21.4,21.8,22.2), mk))   # atteso 21.6 (l'originale dava 11.2)
# dbh_robusto(lapply(c(1,21.8,22.0,22.2), mk))              # atteso 22.0
# dbh_robusto(lapply(c(1,1,1), mk))                         # atteso 1 (ripiego)
# v5 sull'8 a seed di default: prima 16.1 (42.0%), atteso circa 25 cm (errore 7.8-11.7%)


# --- FIX 5: distances_between_lines() - versione fragile da helpers.R ------
# helpers.R definisce una versione NON vettorizzata (assume una singola
# linea, non una matrice di molte linee) che restituisce DistOnLines come
# LISTA (sc, tc) invece di un semplice vettore numerico. parent_cylinder.R
# (funzione centrale, collega ogni segmento alla struttura genitore) fa
# "DistOnLines <= cylinder$length[pc]" aspettandosi un vettore - con la
# lista, questo confronto e' silenziosamente sbagliato (oltre al warning
# di riciclo vettoriale osservato sistematicamente nei log).

distances_between_lines <- function(P1, d1, P2mat, d2mat) {
  if (is.null(dim(P2mat))) P2mat <- matrix(P2mat, nrow = 1)
  if (is.null(dim(d2mat))) d2mat <- matrix(d2mat, nrow = 1)
  n <- nrow(P2mat)
  d1m <- matrix(d1, n, 3, byrow = TRUE)
  P1m <- matrix(P1, n, 3, byrow = TRUE)
  w0 <- P1m - P2mat
  a <- rowSums(d1m * d1m)
  b <- rowSums(d1m * d2mat)
  cc <- rowSums(d2mat * d2mat)
  d <- rowSums(d1m * w0)
  e <- rowSums(d2mat * w0)
  denom <- a * cc - b * b
  denom[abs(denom) < 1e-12] <- 1e-12
  t1 <- (b * e - cc * d) / denom
  t2 <- (a * e - b * d) / denom
  closest1 <- P1m + t1 * d1m
  closest2 <- P2mat + t2 * d2mat
  Dist <- sqrt(rowSums((closest1 - closest2)^2))
  list(Dist = Dist, DistOnLine1 = t1, DistOnLines = t2)
}

# --- FIX 6: create_prioritized_permutation() - bug di sample(x) con --------
# x scalare/lunghezza 1. In R, sample(x) con x un singolo numero N non lo
# tratta come "permuta il vettore {N}" ma come "campiona da 1:N" - un
# comportamento speciale documentato ma facile da dimenticare. Se un
# qualsiasi gruppo (small/medium/large) ha esattamente 1 elemento, la
# funzione produce un vettore di lunghezza completamente sbagliata,
# corrompendo silenziosamente l'ordine di elaborazione dei cover set.
# Verificato sperimentalmente: sample(c(57)) restituisce 57 numeri casuali
# da 1 a 57, non {57}. Fix: usare l'idioma sicuro x[sample.int(length(x))].

create_prioritized_permutation <- function(RelSize) {
  np <- length(RelSize)
  ind <- 1:np
  small <- which(RelSize <= 32)
  medium <- which(RelSize > 32 & RelSize <= 128)
  large <- which(RelSize > 128)
  RandPerm <- integer(np)
  if (length(small) > 0) {
    RandPerm[1:length(small)] <- small[sample.int(length(small))]
  }
  if (length(medium) > 0) {
    start <- length(small) + 1
    RandPerm[start:(start + length(medium) - 1)] <- medium[sample.int(length(medium))]
  }
  if (length(large) > 0) {
    start <- length(small) + length(medium) + 1
    RandPerm[start:np] <- large[sample.int(length(large))]
  }
  return(RandPerm)
}

# ============================================================
# SOSTITUZIONE DEFINITIVA di dbh_robusto() nel modulo che la definisce
# Versione v3: fix base (esclude modelli con DBH < 5 cm) + guard per-modello
#   (se DBHcyl < 0.8 * DBHqsm e DBHqsm >= 5 cm, usa DBHqsm per quel modello).
#
# Decisione presa il 29/09 dopo valutazione multi-seed (16 alberi x 5 offset di seed):
#   - nessun albero peggiora (confronto min/max su 5 offset, tutti e 16)
#   - guard attiva su 2 alberi (25, 3): 25 migliora (mediana 19.6 -> 14.3), 3 quasi invariato (13.9 -> 13.8)
#   - media delle mediane per albero: 11.62 -> 11.29; mediana delle mediane: 9.76 -> 9.76 (invariata)
#   - mediana per offset: 3 offset su 5 migliorano o pari, 2 peggiorano di pochi decimi (dentro il rumore da seed)
# Base statistica stretta (1 solo guadagno chiaro, il 25) ma l'effetto e' motivato dal meccanismo del
#   codice (letto e capito in sessione, non trovato frugando nei risultati) e non produce danni misurati.
#
# Trovare il file con: grep "dbh_robusto <- function" in R/; sostituire nella definizione caricata
# per ultima (nel bootstrap: correzioni_sessione_oggi.R). Dopo la sostituzione non caricare piu'
# fix_dbh_robusto.R ne' dbh_robusto_v3_guard_qsm.R.
# ============================================================

dbh_robusto <- function(QSM_list, min_height_frac = 0.5, mad_soglia = 2.5, min_cilindri = 3,
                        dbh_min_valido = 5, soglia_guard = 0.8) {
  if (length(QSM_list) == 0) return(NA)
  
  # FIX 1 (guard per-modello): se il refit (DBHcyl) e' molto piu' piccolo della lettura della
  # catena di cilindri a 1.3m (DBHqsm) - sintomo di un fit del cerchio finale fallito su un
  # tronco comunque ben segmentato - si usa DBHqsm per quel modello, prima dell'aggregazione.
  QSM_list <- lapply(QSM_list, function(q) {
    if (is.null(q) || is.null(q$treedata)) return(q)
    dcyl <- q$treedata$DBHcyl; dqsm <- q$treedata$DBHqsm
    if (is.null(dcyl) || length(dcyl) == 0 || is.null(dqsm) || length(dqsm) == 0) return(q)
    if (!is.finite(dcyl[1]) || !is.finite(dqsm[1])) return(q)
    if (dqsm[1] * 100 >= dbh_min_valido && dcyl[1] < soglia_guard * dqsm[1]) {
      q$treedata$DBHcyl <- dqsm[1]
    }
    q
  })
  
  # FIX 2 (esclusione dei fit degeneri): una chiamata a treeqsm() produce fino a 8 modelli
  # (griglia di parametri); in alcuni il DBH e' esattamente 1.0 (valore di ripiego di un fit
  # degenere). Aggregarlo come se fosse una misura produce letture a meta' strada o azzerate
  # quando questi modelli sono la meta' o piu' del totale. Si escludono prima dell'aggregazione.
  ok_valido <- vapply(QSM_list, function(q) {
    d <- q$treedata$DBHcyl
    !is.null(d) && length(d) > 0 && is.finite(d[1]) && d[1] * 100 >= dbh_min_valido
  }, logical(1))
  if (any(ok_valido)) QSM_list <- QSM_list[ok_valido]
  
  # --- Aggregazione originale (mediana con filtro a 2.5 MAD) ---
  non_null <- sapply(QSM_list, function(q) !is.null(q) && !is.null(q$treedata))
  if (!any(non_null)) return(NA)
  heights <- sapply(QSM_list[non_null], function(q) {
    h <- q$treedata$TreeHeight
    if (is.null(h) || length(h) == 0 || !is.finite(h)) NA_real_ else h
  })
  med_h <- median(heights, na.rm = TRUE)
  ok_height <- !is.na(heights) & (is.na(med_h) | heights >= min_height_frac * med_h)
  n_cyl <- sapply(QSM_list[non_null], function(q) length(q$cylinder$radius))
  ok_cyl <- n_cyl >= min_cilindri
  dbh_grezzo <- sapply(QSM_list[non_null], function(q) {
    d <- q$treedata$DBHcyl
    if (is.null(d) || length(d) == 0) NA_real_ else d
  })
  ok_dbh_finito <- is.finite(dbh_grezzo)
  ok_tutti <- ok_height & ok_cyl & ok_dbh_finito
  candidates_idx <- which(non_null)[ok_tutti]
  if (length(candidates_idx) == 0) return(NA)
  dbh_vals <- sapply(QSM_list[candidates_idx], function(q) q$treedata$DBHcyl * 100)
  med_dbh <- median(dbh_vals, na.rm = TRUE)
  mad_dbh <- mad(dbh_vals, na.rm = TRUE)
  if (is.na(mad_dbh) || mad_dbh < 1e-6) return(med_dbh)
  ok_dbh <- abs(dbh_vals - med_dbh) <= mad_soglia * mad_dbh
  if (sum(ok_dbh) == 0) return(med_dbh)
  median(dbh_vals[ok_dbh], na.rm = TRUE)
}

# --- Verifiche dopo la sostituzione ---
# mk <- function(dcyl, dqsm, h=10) list(treedata=list(DBHcyl=dcyl/100, DBHqsm=dqsm/100, TreeHeight=h), cylinder=list(radius=rep(0.1,50)))
# dbh_robusto(list(mk(14.4, 24.1)))                                    # atteso 24.1 (guard scatta)
# dbh_robusto(lapply(c(21.4,21.4,21.8,22.2), function(x) mk(x,x)))     # atteso 21.6-22.0 (invariato)
# dbh_robusto(list(mk(1,3)))                                           # atteso 1 (DBHqsm<5, guard non scatta)

# --- PIPELINE FINALE (Day 22): ricentraggio + controllo di plausibilita' ---
# Versione che ha prodotto risultati_v20: mediana 20.5%, media 27.18%
# (vedi riepilogo, sezione "Day 22" per dettagli e tabella completa)

pipeline_integrata_v11 <- function(P_grezza_input, inputs, inputs_vecchi, n_modelli = 5,
                                   soglia_calo = 0.07, soglia_residui = 0.075,
                                   max_iter = 6, raggio_max_plausibile = 80,
                                   soglia_spostamento_centro = 0.7, soglia_ricentraggio = 0.5,
                                   soglia_restringimento = 0.5) {
  
  P_grezza <- P_grezza_input
  disp <- rileva_dispersione_anomala(P_grezza)
  P_pulito <- P_grezza[P_grezza[,3] >= disp$h_inizio_tronco_pulito, , drop = FALSE]
  q <- treeqsm(P_pulito, inputs)
  cyl <- modello_non_degenere(q)$cylinder
  
  is_tronco_r <- cyl$branch == 1
  x_tronco <- median(cyl$start[is_tronco_r, 1]); y_tronco <- median(cyl$start[is_tronco_r, 2])
  spostamento_iniziale <- sqrt(x_tronco^2 + y_tronco^2)
  ricentrato <- FALSE
  if (spostamento_iniziale > soglia_ricentraggio) {
    P_grezza[,1] <- P_grezza[,1] - x_tronco; P_grezza[,2] <- P_grezza[,2] - y_tronco
    P_pulito[,1] <- P_pulito[,1] - x_tronco; P_pulito[,2] <- P_pulito[,2] - y_tronco
    ricentrato <- TRUE
    q <- treeqsm(P_pulito, inputs)
    cyl <- modello_non_degenere(q)$cylinder
  }
  
  dbh_standard <- dbh_robusto(q)
  if (is.na(dbh_standard)) dbh_standard <- q[[1]]$treedata$DBHcyl * 100
  
  diag <- diagnostica_rumore_interno_diretta(P_grezza)
  e_catena <- diagnostica_catena_discontinua(cyl)
  
  ref_base <- robust_trunk_reference(cyl, cyl$branch==1, seq_along(cyl$radius))
  dist_centro_iniziale <- sqrt(ref_base$start[1]^2 + ref_base$start[2]^2)
  usa_fascia_bassa <- dist_centro_iniziale > soglia_spostamento_centro
  
  applica_knn <- e_catena || (!is.null(diag$esito) && diag$esito == "RUMORE_INTERNO_DIFFUSO" && !usa_fascia_bassa)
  
  if (applica_knn) {
    keep_knn <- tryCatch(filtro_knn_fedele(P_grezza), error=function(e) rep(TRUE, nrow(P_grezza)))
    P_grezza_knn <- P_grezza[keep_knn, , drop = FALSE]
    disp2 <- tryCatch(rileva_dispersione_anomala(P_grezza_knn), error=function(e) NULL)
    if (!is.null(disp2)) {
      P_pulito2 <- P_grezza_knn[P_grezza_knn[,3] >= disp2$h_inizio_tronco_pulito, , drop = FALSE]
      q2 <- tryCatch(treeqsm(P_pulito2, inputs), error=function(e) NULL)
      if (!is.null(q2)) {
        P_grezza <- P_grezza_knn; P_pulito <- P_pulito2
        cyl <- modello_non_degenere(q2)$cylinder
        dbh_standard2 <- dbh_robusto(q2)
        dbh_standard <- if (is.na(dbh_standard2)) q2[[1]]$treedata$DBHcyl * 100 else dbh_standard2
      }
    }
  }
  
  risultato <- list(DBHcyl_finale = dbh_standard, DBHcyl_standard = dbh_standard,
                    diagnosi = diag$esito, metodo = "standard", ricentrato = ricentrato)
  
  if (!is.null(diag$esito) && diag$esito == "RUMORE_INTERNO_DIFFUSO") {
    dbh_corretti <- c()
    for (m in seq_len(n_modelli)) {
      inputs_m <- inputs; inputs_m$model <- m
      q_m <- tryCatch(treeqsm(P_pulito, inputs_m), error = function(e) NULL)
      if (is.null(q_m)) next
      inputs_m_flat <- appiattisci_inputs(inputs_m)
      cyl_m <- modello_non_degenere(q_m)$cylinder
      ref_m <- robust_trunk_reference(cyl_m, cyl_m$branch==1, seq_along(cyl_m$radius))
      Z_start_m <- ref_m$start[3]
      cover1_m <- cover_sets(P_pulito, inputs_m_flat)
      ts1_m <- tree_sets(P_pulito, cover1_m, inputs_m_flat)
      seg1_m <- correct_segments(P_pulito, ts1_m$cover, segments(ts1_m$cover, ts1_m$Base, ts1_m$Forb), inputs_m_flat,0,1,1)
      RS_m <- relative_size(P_pulito, ts1_m$cover, seg1_m)
      cover2_m <- cover_sets(P_pulito, inputs_m_flat, RS_m)
      ts2_m <- tree_sets(P_pulito, cover2_m, inputs_m_flat, seg1_m)
      seg2_m <- correct_segments(P_pulito, ts2_m$cover, segments(ts2_m$cover, ts2_m$Base, ts2_m$Forb), inputs_m_flat,1,1,0)
      trunk_idx_m <- unique(unlist(ts2_m$cover$ball[unlist(seg2_m$segments[[1]])]))
      trunk_pts_m <- P_pulito[trunk_idx_m, , drop = FALSE]
      V <- sweep(trunk_pts_m,2,ref_m$start); hh_m <- as.vector(V %*% ref_m$axis)
      if (Z_start_m > 1.35) { fascia_ass <- c(Z_start_m+0.02, Z_start_m+0.42) }
      else if (usa_fascia_bassa) { fascia_ass <- c(1.2, 1.4) } else { fascia_ass <- c(1.7, 2.1) }
      fascia_hh <- fascia_ass - Z_start_m
      if (fascia_hh[2] <= 0.05) next
      Tpts <- trunk_pts_m[hh_m > fascia_hh[1] & hh_m < fascia_hh[2], , drop = FALSE]
      if (nrow(Tpts) < 50) next
      dl_init <- distances_to_line(Tpts, ref_m$axis, ref_m$start)
      km <- tryCatch(kmeans(dl_init$d, centers = 2), error = function(e) NULL)
      if (is.null(km)) next
      mask_ext <- km$cluster == which.max(km$centers)
      dbh_corretti <- c(dbh_corretti, mean(dl_init$d[mask_ext])*2*100)
    }
    if (length(dbh_corretti) > 0) { risultato$DBHcyl_finale <- median(dbh_corretti); risultato$metodo <- "multimodello_centroguidato" }
    return(risultato)
  }
  
  if (!is.null(diag$esito) && diag$esito == "contaminazione_localizzata") {
    P_taglio05 <- P_grezza[P_grezza[,3] >= 0.5, , drop = FALSE]
    q_t05 <- tryCatch(treeqsm(P_taglio05, inputs_vecchi), error = function(e) NULL)
    if (!is.null(q_t05)) {
      dbh_t05 <- dbh_robusto(q_t05)
      if (!is.na(dbh_t05)) { risultato$DBHcyl_finale <- dbh_t05; risultato$metodo <- "taglio05_vecchi_param" }
      else { risultato$metodo <- "taglio05_fallito_tutti_degeneri" }
    }
    due_fusti <- tryCatch(traccia_due_alberi_da_altezza(P_grezza), error = function(e) NULL)
    if (!is.null(due_fusti)) {
      A_t <- due_fusti$albero_A[due_fusti$albero_A[,3] >= 0.5, , drop = FALSE]
      B_t <- due_fusti$albero_B[due_fusti$albero_B[,3] >= 0.5, , drop = FALSE]
      qA <- tryCatch(treeqsm(A_t, inputs_vecchi), error = function(e) NULL)
      qB <- tryCatch(treeqsm(B_t, inputs_vecchi), error = function(e) NULL)
      if (!is.null(qA) && !is.null(qB)) {
        dA <- dbh_robusto(qA); dB <- dbh_robusto(qB)
        if (!is.na(dA) && !is.na(dB) && min(dA,dB)/max(dA,dB) > 0.75) {
          risultato$DBHcyl_finale <- sqrt(dA^2+dB^2); risultato$metodo <- "diametro_equivalente"
        }
      }
    }
    return(risultato)
  }
  
  if (e_catena) {
    is_tronco0 <- cyl$branch == 1
    r_c0 <- cyl$radius[is_tronco0] * 100
    km_cat <- tryCatch(kmeans(r_c0, centers = 2), error = function(e) NULL)
    if (!is.null(km_cat)) {
      cluster_alto0 <- which.max(km_cat$centers)
      risultato$DBHcyl_finale <- mean(r_c0[km_cat$cluster == cluster_alto0]) * 2
      risultato$diagnosi <- "catena_discontinua"; risultato$metodo <- "clustering_raggio_catena"
      return(risultato)
    }
  }
  
  inputs_flat <- appiattisci_inputs(inputs)
  ref <- robust_trunk_reference(cyl, cyl$branch==1, seq_along(cyl$radius))
  Z_start <- ref$start[3]
  cover1 <- cover_sets(P_pulito, inputs_flat)
  ts1 <- tree_sets(P_pulito, cover1, inputs_flat)
  seg1 <- correct_segments(P_pulito, ts1$cover, segments(ts1$cover, ts1$Base, ts1$Forb), inputs_flat,0,1,1)
  RS <- relative_size(P_pulito, ts1$cover, seg1)
  cover2 <- cover_sets(P_pulito, inputs_flat, RS)
  ts2 <- tree_sets(P_pulito, cover2, inputs_flat, seg1)
  seg2 <- correct_segments(P_pulito, ts2$cover, segments(ts2$cover, ts2$Base, ts2$Forb), inputs_flat,1,1,0)
  trunk_idx <- unique(unlist(ts2$cover$ball[unlist(seg2$segments[[1]])]))
  trunk_pts <- P_pulito[trunk_idx, , drop = FALSE]
  V <- sweep(trunk_pts, 2, ref$start); hh <- as.vector(V %*% ref$axis)
  fasce_assolute <- list(c(1.2,1.4), c(1.1,1.6), c(0.9,1.4))
  if (Z_start > 1.35) fasce_assolute <- c(fasce_assolute, list(c(Z_start+0.02, Z_start+0.42)))
  migliore_sd <- Inf
  for (fascia_ass in fasce_assolute) {
    fascia_hh <- fascia_ass - Z_start
    if (fascia_hh[2] <= 0.05) next
    mask_zona <- hh > fascia_hh[1] & hh < fascia_hh[2]
    pts_zona <- trunk_pts[mask_zona, , drop = FALSE]
    if (nrow(pts_zona) < 30) next
    d_zona <- distances_to_line(pts_zona, ref$axis, ref$start)$d * 100
    mediana_grezza <- median(d_zona)
    if (mediana_grezza > raggio_max_plausibile) next
    idx_correnti <- seq_along(d_zona); rapporto_prec <- NA; idx_migliori <- idx_correnti
    for (iter in 1:max_iter) {
      km <- tryCatch(kmeans(d_zona[idx_correnti], centers = 2), error = function(e) NULL)
      if (is.null(km)) break
      centri <- sort(km$centers); rapporto <- centri[2]/centri[1]
      cluster_int <- which.min(km$centers); idx_nuovi <- idx_correnti[km$cluster == cluster_int]
      calo <- if(!is.na(rapporto_prec)) (rapporto_prec - rapporto)/rapporto_prec else NA
      if (!is.na(calo) && calo > soglia_calo) break
      if (length(idx_nuovi) < 20) break
      idx_migliori <- idx_nuovi; idx_correnti <- idx_nuovi; rapporto_prec <- rapporto
    }
    d_migliori <- d_zona[idx_migliori]
    if (mean(d_migliori) < soglia_restringimento * mediana_grezza) next
    pts_finali <- pts_zona[idx_migliori, , drop = FALSE]
    cyl0 <- list(radius = mean(d_migliori)/100, start = ref$start, axis = ref$axis)
    cyl_fit <- tryCatch(least_squares_cylinder(pts_finali, cyl0), error = function(e) NULL)
    if (is.null(cyl_fit) || !isTRUE(cyl_fit$conv)) next
    d_residui <- distances_to_line(pts_finali, cyl_fit$axis, cyl_fit$start)$d * 100
    sd_relativo <- sd(d_residui - cyl_fit$radius*100) / (cyl_fit$radius*100)
    if (sd_relativo < migliore_sd && sd_relativo <= soglia_residui) {
      migliore_sd <- sd_relativo
      risultato$DBHcyl_finale <- mean(d_migliori) * 2
      risultato$metodo <- "clustering_universale"
      risultato$sd_relativo <- round(sd_relativo,3)
      risultato$fascia_usata <- paste(round(fascia_ass,2), collapse="-")
    }
  }
  risultato
}

# --- NUOVO STRUMENTO (Day 22): criterio multifascia per il raggio di --------
# ri-estrazione sicuro dal LAS, per alberi con point cloud pre-estratta
# gravemente povera di punti (categoria "contaminazione_localizzata").
#
# Idea (suggerita dall'utente): invece di scegliere un raggio di taglio
# fisso o cercarlo per tentativi, testare un clustering/gap radiale su
# PIU' fasce di altezza contemporaneamente e usare la mediana dei raggi-gap
# trovati sulle fasce con segnale sufficiente. Un gap ampio (>0.5m di norma)
# in una fascia indica il vero confine tra il tronco (visto solo in parte,
# spesso per occlusione angolare) e un oggetto esterno (albero vicino).
#
# REQUISITO IMPORTANTE, scoperto empiricamente: la point cloud sorgente
# (da ri-estrarre dal LAS) deve avere un raggio di estrazione GENEROSO
# (5m verificato, 3m insufficiente) per dare al criterio abbastanza fasce
# con segnale valido. Con pochi punti, spesso emerge una sola fascia
# valida (o nessuna), e il risultato diventa fragile.
#
# Verificato su due casi (entrambi con raggio di estrazione LAS = 5m):
# - Albero 20 (treeID 20): pipeline attuale 96.1% di errore -> 0.6%
#   (6 fasce concordanti nel test con raggio estrazione 3m, risultato gia'
#   ottimo; non riverificato a 5m dato il risultato gia' eccellente)
# - Albero 25/MATLAB (treeID 38): pipeline attuale 41.4% -> 18.7% con
#   estrazione a 3m (1 sola fascia valida, segnale debole) -> 7.5% con
#   estrazione a 5m (6 fasce concordanti, segnale robusto)
#
# NON ANCORA INTEGRATO nella pipeline di produzione: verificato solo su
# 2 casi, entrambi nella categoria "contaminazione_localizzata". Il flusso
# completo per usarlo: 1) leggere il LAS con raggio 5m attorno al centro
# nominale dell'albero, 2) lidR::classify_ground(csf()) + normalize_height
# (nessun punto di terreno pre-classificato in questi LAS), 3) centrare
# sul centro nominale, 4) chiamare questa funzione per trovare il raggio
# sicuro, 5) tagliare la point cloud a quel raggio e rifare il fit standard.
#
# @param P Point cloud (matrice n x 3, gia' centrata sul centro nominale
#   e con Z normalizzata sul terreno)
# @param h_max Altezza massima da considerare (default: min(8, max altezza))
# @param larghezza_fascia Ampiezza di ciascuna fascia testata, in metri
# @param ampiezza_min_gap Ampiezza minima del gap perche' una fascia sia
#   considerata valida (scarta fasce con solo rumore locale)
# @param raggio_min, raggio_max Range di raggi considerati per il gap
# @param margine_sicurezza Fattore moltiplicativo di sicurezza (<1) sul
#   raggio trovato, per non tagliare esattamente al bordo del gap
# @return Lista con raggio_sicuro (NA se nessuna fascia valida),
#   n_fasce_valide (indicatore di affidabilita' del risultato: piu' alto
#   = piu' robusto, verificato empiricamente sui due casi sopra),
#   raggi_per_fascia (i singoli valori usati per la mediana)

trova_raggio_sicuro_multifascia <- function(P, h_max = NULL, larghezza_fascia = 1.0, 
                                             ampiezza_min_gap = 0.5, raggio_min = 0.05, raggio_max = 5,
                                             margine_sicurezza = 0.95) {
  if (is.null(h_max)) h_max <- min(8, max(P[,3]))
  fasce <- seq(0.9, h_max - larghezza_fascia, by = larghezza_fascia)
  raggi_trovati <- c()
  for (h0 in fasce) {
    mask_z <- P[,3] >= h0 & P[,3] < h0 + larghezza_fascia
    pts_z <- P[mask_z, 1:2, drop=FALSE]
    if (nrow(pts_z) < 30) next
    d <- sqrt(pts_z[,1]^2 + pts_z[,2]^2)
    d <- sort(d[d >= raggio_min & d <= raggio_max])
    if (length(d) < 30) next
    gaps <- diff(d)
    idx_max <- which.max(gaps)
    if (gaps[idx_max] >= ampiezza_min_gap) {
      raggi_trovati <- c(raggi_trovati, d[idx_max])
    }
  }
  if (length(raggi_trovati) == 0) return(list(raggio_sicuro = NA, n_fasce_valide = 0))
  list(raggio_sicuro = median(raggi_trovati) * margine_sicurezza, 
       n_fasce_valide = length(raggi_trovati), 
       raggi_per_fascia = raggi_trovati)
}

# --- NUOVO STRUMENTO (Day 22, seconda parte): funzione unificata per la ---
# ri-estrazione dal LAS + criterio multifascia, pronta per la pipeline.
# Verificato su 2 casi: albero 20 (~1.3% errore) e albero 38/25 MATLAB
# (~7.4% errore), entrambi con lo stesso identico processo automatico
# (nessun parametro scelto caso per caso).
#
# @param cx, cy Coordinate nominali dell'albero nel sistema del LAS
# @param file_las Vettore di path ai file LAS disponibili (la funzione
#   sceglie automaticamente quale copre le coordinate date)
# @param raggio_estrazione Raggio (m) per la lettura iniziale dal LAS
#   (default 5m, verificato necessario - 3m insufficiente in un caso)
# @param max_punti_fit Punti massimi passati a treeqsm() (campionamento
#   casuale se superati, per velocita')
# @return Lista con successo (TRUE/FALSE), e se TRUE: P (la point cloud
#   pronta per il fit), raggio_usato, n_fasce_valide (indicatore di
#   affidabilita'), n_punti, las_usato. Se FALSE: motivo del fallimento.

file_las_disponibili <- c(
  "F:/EB/tartarino/sw python/TreeQSM_ProjectR/Data/Area_2_LAS_11.las",
  "F:/EB/tartarino/sw python/TreeQSM_ProjectR/Data/Area_2_LAS_15.las"
)

arricchisci_da_las <- function(cx, cy, file_las = file_las_disponibili, raggio_estrazione = 5,
                                max_punti_fit = 15000, seed = 42) {
  las_scelto <- NULL
  for (f in file_las) {
    hdr <- tryCatch(lidR::readLASheader(f), error = function(e) NULL)
    if (is.null(hdr)) next
    if (cx >= hdr@PHB[["Min X"]] && cx <= hdr@PHB[["Max X"]] &&
        cy >= hdr@PHB[["Min Y"]] && cy <= hdr@PHB[["Max Y"]]) {
      las_scelto <- f; break
    }
  }
  if (is.null(las_scelto)) return(list(successo = FALSE, motivo = "nessun LAS copre queste coordinate"))
  las <- tryCatch(lidR::readLAS(las_scelto, 
                    filter = sprintf("-drop_x_below %f -drop_x_above %f -drop_y_below %f -drop_y_above %f",
                                      cx - raggio_estrazione, cx + raggio_estrazione,
                                      cy - raggio_estrazione, cy + raggio_estrazione)),
                  error = function(e) NULL)
  if (is.null(las) || nrow(las@data) < 100) return(list(successo = FALSE, motivo = "lettura LAS fallita o troppo pochi punti"))
  las_ground <- tryCatch(lidR::classify_ground(las, lidR::csf()), error = function(e) NULL)
  if (is.null(las_ground)) return(list(successo = FALSE, motivo = "classificazione terreno fallita"))
  las_norm <- tryCatch(lidR::normalize_height(las_ground, lidR::knnidw(k = 10)), error = function(e) NULL)
  if (is.null(las_norm)) return(list(successo = FALSE, motivo = "normalizzazione altezza fallita"))
  P <- as.matrix(data.frame(X = las_norm@data$X - cx, Y = las_norm@data$Y - cy, Z = las_norm@data$Z))
  P <- P[!is.na(P[,3]) & P[,3] >= 0, , drop = FALSE]
  if (nrow(P) < 100) return(list(successo = FALSE, motivo = "troppo pochi punti dopo pulizia"))
  ris_gap <- trova_raggio_sicuro_multifascia(P)
  if (is.na(ris_gap$raggio_sicuro)) return(list(successo = FALSE, motivo = "nessuna fascia valida trovata"))
  mask_finale <- sqrt(P[,1]^2 + P[,2]^2) < ris_gap$raggio_sicuro
  P_finale <- P[mask_finale, , drop = FALSE]
  if (nrow(P_finale) < 100) return(list(successo = FALSE, motivo = "troppo pochi punti dopo il taglio finale"))
  set.seed(seed)
  if (nrow(P_finale) > max_punti_fit) {
    P_finale <- P_finale[sample(nrow(P_finale), max_punti_fit), , drop = FALSE]
  }
  list(successo = TRUE, P = P_finale, raggio_usato = ris_gap$raggio_sicuro, 
       n_fasce_valide = ris_gap$n_fasce_valide, n_punti = nrow(P_finale), las_usato = basename(las_scelto))
}

# --- AGGIORNAMENTO (Day 22, terza parte): criterio di affidabilita' -------
# rivisto. Il coefficiente di variazione tra le fasce (usato in un primo
# tentativo) non discriminava bene: l'albero 20 (buon caso, 1.6% errore
# nel test finale) aveva CV=0.348, l'albero 17 (cattivo caso, 96.4% errore)
# CV=0.379 - margine troppo stretto per essere affidabile su soli 3 casi.
#
# Criterio sostituito con un limite fisico assoluto sul raggio trovato:
# un vero tronco, anche gravemente occluso, non dovrebbe richiedere un
# raggio di taglio superiore a ~2m. Un raggio piu' grande indica che il
# criterio ha probabilmente gia' agganciato l'oggetto esterno invece del
# vero tronco. Verificato: distingue correttamente tutti e 3 i casi noti
# (20 e 38: raggio ~1.1m, affidabile; 17: raggio 3.85m, non affidabile).
#
# Nuovo campo "affidabile" (TRUE/FALSE): combina n_fasce_valide >= 3 E
# raggio_sicuro <= raggio_max_plausibile_tronco. Il chiamante deve
# controllare questo flag prima di usare il risultato, non solo
# raggio_sicuro o n_fasce_valide singolarmente.

trova_raggio_sicuro_multifascia <- function(P, h_max = NULL, larghezza_fascia = 1.0, 
                                             ampiezza_min_gap = 0.5, raggio_min = 0.05, raggio_max = 5,
                                             margine_sicurezza = 0.95, raggio_max_plausibile_tronco = 2.0,
                                             n_fasce_min = 3) {
  if (is.null(h_max)) h_max <- min(8, max(P[,3]))
  fasce <- seq(0.9, h_max - larghezza_fascia, by = larghezza_fascia)
  raggi_trovati <- c()
  for (h0 in fasce) {
    mask_z <- P[,3] >= h0 & P[,3] < h0 + larghezza_fascia
    pts_z <- P[mask_z, 1:2, drop=FALSE]
    if (nrow(pts_z) < 30) next
    d <- sqrt(pts_z[,1]^2 + pts_z[,2]^2)
    d <- sort(d[d >= raggio_min & d <= raggio_max])
    if (length(d) < 30) next
    gaps <- diff(d)
    idx_max <- which.max(gaps)
    if (gaps[idx_max] >= ampiezza_min_gap) {
      raggi_trovati <- c(raggi_trovati, d[idx_max])
    }
  }
  if (length(raggi_trovati) == 0) return(list(raggio_sicuro = NA, n_fasce_valide = 0, affidabile = FALSE))
  raggio_finale <- median(raggi_trovati) * margine_sicurezza
  affidabile <- length(raggi_trovati) >= n_fasce_min && raggio_finale <= raggio_max_plausibile_tronco
  list(raggio_sicuro = raggio_finale, n_fasce_valide = length(raggi_trovati), 
       affidabile = affidabile, raggi_per_fascia = raggi_trovati)
}

# --- FIX 8: arricchisci_da_las() non propagava "affidabile" da ------------
# trova_raggio_sicuro_multifascia(), e usava il nome "raggio_usato" invece
# di "raggio_sicuro" (incoerente con il resto). La pipeline cercava
# ris_las$affidabile (sempre NULL/FALSE) e ris_las$raggio_sicuro (sempre
# NULL) - risultato: MAI usava il ramo ad alta affidabilita' (fallback
# sempre attivato), e in piu' un errore vero e proprio su round(NULL,2)
# quando si costruiva il messaggio di avviso. Corretto qui sotto: nomi
# coerenti, "affidabile" propagato correttamente.

arricchisci_da_las <- function(cx, cy, file_las = file_las_disponibili, raggio_estrazione = 5,
                                max_punti_fit = 15000, seed = 42) {
  las_scelto <- NULL
  for (f in file_las) {
    hdr <- tryCatch(lidR::readLASheader(f), error = function(e) NULL)
    if (is.null(hdr)) next
    if (cx >= hdr@PHB[["Min X"]] && cx <= hdr@PHB[["Max X"]] &&
        cy >= hdr@PHB[["Min Y"]] && cy <= hdr@PHB[["Max Y"]]) {
      las_scelto <- f; break
    }
  }
  if (is.null(las_scelto)) return(list(successo = FALSE, motivo = "nessun LAS copre queste coordinate"))
  las <- tryCatch(lidR::readLAS(las_scelto, 
                    filter = sprintf("-drop_x_below %f -drop_x_above %f -drop_y_below %f -drop_y_above %f",
                                      cx - raggio_estrazione, cx + raggio_estrazione,
                                      cy - raggio_estrazione, cy + raggio_estrazione)),
                  error = function(e) NULL)
  if (is.null(las) || nrow(las@data) < 100) return(list(successo = FALSE, motivo = "lettura LAS fallita o troppo pochi punti"))
  las_ground <- tryCatch(lidR::classify_ground(las, lidR::csf()), error = function(e) NULL)
  if (is.null(las_ground)) return(list(successo = FALSE, motivo = "classificazione terreno fallita"))
  las_norm <- tryCatch(lidR::normalize_height(las_ground, lidR::knnidw(k = 10)), error = function(e) NULL)
  if (is.null(las_norm)) return(list(successo = FALSE, motivo = "normalizzazione altezza fallita"))
  P <- as.matrix(data.frame(X = las_norm@data$X - cx, Y = las_norm@data$Y - cy, Z = las_norm@data$Z))
  P <- P[!is.na(P[,3]) & P[,3] >= 0, , drop = FALSE]
  if (nrow(P) < 100) return(list(successo = FALSE, motivo = "troppo pochi punti dopo pulizia"))
  ris_gap <- trova_raggio_sicuro_multifascia(P)
  if (is.na(ris_gap$raggio_sicuro)) return(list(successo = FALSE, motivo = "nessuna fascia valida trovata"))
  mask_finale <- sqrt(P[,1]^2 + P[,2]^2) < ris_gap$raggio_sicuro
  P_finale <- P[mask_finale, , drop = FALSE]
  if (nrow(P_finale) < 100) return(list(successo = FALSE, motivo = "troppo pochi punti dopo il taglio finale"))
  set.seed(seed)
  if (nrow(P_finale) > max_punti_fit) {
    P_finale <- P_finale[sample(nrow(P_finale), max_punti_fit), , drop = FALSE]
  }
  list(successo = TRUE, P = P_finale, raggio_sicuro = ris_gap$raggio_sicuro, 
       n_fasce_valide = ris_gap$n_fasce_valide, affidabile = ris_gap$affidabile,
       n_punti = nrow(P_finale), las_usato = basename(las_scelto))
}

# --- PIPELINE FINALE DEFINITIVA (Day 22, chiusura sessione) ----------------
# Versione che ha prodotto risultati_v23: mediana 8.8%, media 18.23%
# sui 15 alberi con risultato (il 42 fallisce sempre per pochi punti).
# Sostituisce la versione precedente in questo stesso file (Day 22,
# "PIPELINE FINALE": ricentraggio + controllo di plausibilita') aggiungendo
# l'integrazione con arricchisci_da_las() nel ramo contaminazione_localizzata.
# Poiche' i file vengono caricati in ordine e l'ultima definizione vince,
# questa e' la versione che verra' effettivamente usata.
#
# NUOVO rispetto alla precedente: parametri cx_nominale, cy_nominale
# (opzionali - se non forniti, il ramo contaminazione_localizzata usa
# solo il fallback originale, come prima). Se forniti, tenta prima il
# ripopolamento dal LAS; se il risultato non e' affidabile (vedi
# arricchisci_da_las()), fa fallback al metodo originale E imposta
# risultato$avviso con un messaggio esplicito di incertezza residua.

pipeline_integrata_v11 <- function(P_grezza_input, inputs, inputs_vecchi, n_modelli = 5,
                                   soglia_calo = 0.07, soglia_residui = 0.075,
                                   max_iter = 6, raggio_max_plausibile = 80,
                                   soglia_spostamento_centro = 0.7, soglia_ricentraggio = 0.5,
                                   soglia_restringimento = 0.5, cx_nominale = NULL, cy_nominale = NULL) {
  
  P_grezza <- P_grezza_input
  disp <- rileva_dispersione_anomala(P_grezza)
  P_pulito <- P_grezza[P_grezza[,3] >= disp$h_inizio_tronco_pulito, , drop = FALSE]
  q <- treeqsm(P_pulito, inputs)
  cyl <- modello_non_degenere(q)$cylinder
  
  is_tronco_r <- cyl$branch == 1
  x_tronco <- median(cyl$start[is_tronco_r, 1]); y_tronco <- median(cyl$start[is_tronco_r, 2])
  spostamento_iniziale <- sqrt(x_tronco^2 + y_tronco^2)
  ricentrato <- FALSE
  if (spostamento_iniziale > soglia_ricentraggio) {
    P_grezza[,1] <- P_grezza[,1] - x_tronco; P_grezza[,2] <- P_grezza[,2] - y_tronco
    P_pulito[,1] <- P_pulito[,1] - x_tronco; P_pulito[,2] <- P_pulito[,2] - y_tronco
    ricentrato <- TRUE
    q <- treeqsm(P_pulito, inputs)
    cyl <- modello_non_degenere(q)$cylinder
  }
  
  dbh_standard <- dbh_robusto(q)
  if (is.na(dbh_standard)) dbh_standard <- q[[1]]$treedata$DBHcyl * 100
  
  diag <- diagnostica_rumore_interno_diretta(P_grezza)
  e_catena <- diagnostica_catena_discontinua(cyl)
  
  ref_base <- robust_trunk_reference(cyl, cyl$branch==1, seq_along(cyl$radius))
  dist_centro_iniziale <- sqrt(ref_base$start[1]^2 + ref_base$start[2]^2)
  usa_fascia_bassa <- dist_centro_iniziale > soglia_spostamento_centro
  
  applica_knn <- e_catena || (!is.null(diag$esito) && diag$esito == "RUMORE_INTERNO_DIFFUSO" && !usa_fascia_bassa)
  
  if (applica_knn) {
    keep_knn <- tryCatch(filtro_knn_fedele(P_grezza), error=function(e) rep(TRUE, nrow(P_grezza)))
    P_grezza_knn <- P_grezza[keep_knn, , drop = FALSE]
    disp2 <- tryCatch(rileva_dispersione_anomala(P_grezza_knn), error=function(e) NULL)
    if (!is.null(disp2)) {
      P_pulito2 <- P_grezza_knn[P_grezza_knn[,3] >= disp2$h_inizio_tronco_pulito, , drop = FALSE]
      q2 <- tryCatch(treeqsm(P_pulito2, inputs), error=function(e) NULL)
      if (!is.null(q2)) {
        P_grezza <- P_grezza_knn; P_pulito <- P_pulito2
        cyl <- modello_non_degenere(q2)$cylinder
        dbh_standard2 <- dbh_robusto(q2)
        dbh_standard <- if (is.na(dbh_standard2)) q2[[1]]$treedata$DBHcyl * 100 else dbh_standard2
      }
    }
  }
  
  risultato <- list(DBHcyl_finale = dbh_standard, DBHcyl_standard = dbh_standard,
                    diagnosi = diag$esito, metodo = "standard", ricentrato = ricentrato, avviso = NULL)
  
  if (!is.null(diag$esito) && diag$esito == "RUMORE_INTERNO_DIFFUSO") {
    dbh_corretti <- c()
    for (m in seq_len(n_modelli)) {
      inputs_m <- inputs; inputs_m$model <- m
      q_m <- tryCatch(treeqsm(P_pulito, inputs_m), error = function(e) NULL)
      if (is.null(q_m)) next
      inputs_m_flat <- appiattisci_inputs(inputs_m)
      cyl_m <- modello_non_degenere(q_m)$cylinder
      ref_m <- robust_trunk_reference(cyl_m, cyl_m$branch==1, seq_along(cyl_m$radius))
      Z_start_m <- ref_m$start[3]
      cover1_m <- cover_sets(P_pulito, inputs_m_flat)
      ts1_m <- tree_sets(P_pulito, cover1_m, inputs_m_flat)
      seg1_m <- correct_segments(P_pulito, ts1_m$cover, segments(ts1_m$cover, ts1_m$Base, ts1_m$Forb), inputs_m_flat,0,1,1)
      RS_m <- relative_size(P_pulito, ts1_m$cover, seg1_m)
      cover2_m <- cover_sets(P_pulito, inputs_m_flat, RS_m)
      ts2_m <- tree_sets(P_pulito, cover2_m, inputs_m_flat, seg1_m)
      seg2_m <- correct_segments(P_pulito, ts2_m$cover, segments(ts2_m$cover, ts2_m$Base, ts2_m$Forb), inputs_m_flat,1,1,0)
      trunk_idx_m <- unique(unlist(ts2_m$cover$ball[unlist(seg2_m$segments[[1]])]))
      trunk_pts_m <- P_pulito[trunk_idx_m, , drop = FALSE]
      V <- sweep(trunk_pts_m,2,ref_m$start); hh_m <- as.vector(V %*% ref_m$axis)
      if (Z_start_m > 1.35) { fascia_ass <- c(Z_start_m+0.02, Z_start_m+0.42) }
      else if (usa_fascia_bassa) { fascia_ass <- c(1.2, 1.4) } else { fascia_ass <- c(1.7, 2.1) }
      fascia_hh <- fascia_ass - Z_start_m
      if (fascia_hh[2] <= 0.05) next
      Tpts <- trunk_pts_m[hh_m > fascia_hh[1] & hh_m < fascia_hh[2], , drop = FALSE]
      if (nrow(Tpts) < 50) next
      dl_init <- distances_to_line(Tpts, ref_m$axis, ref_m$start)
      km <- tryCatch(kmeans(dl_init$d, centers = 2), error = function(e) NULL)
      if (is.null(km)) next
      mask_ext <- km$cluster == which.max(km$centers)
      dbh_corretti <- c(dbh_corretti, mean(dl_init$d[mask_ext])*2*100)
    }
    if (length(dbh_corretti) > 0) { risultato$DBHcyl_finale <- median(dbh_corretti); risultato$metodo <- "multimodello_centroguidato" }
    return(risultato)
  }
  
  if (!is.null(diag$esito) && diag$esito == "contaminazione_localizzata") {
    ris_las <- NULL
    if (!is.null(cx_nominale) && !is.null(cy_nominale)) {
      ris_las <- tryCatch(arricchisci_da_las(cx_nominale, cy_nominale), error = function(e) list(successo = FALSE, motivo = conditionMessage(e)))
    }
    
    if (!is.null(ris_las) && isTRUE(ris_las$successo) && isTRUE(ris_las$affidabile)) {
      q_las <- tryCatch(treeqsm(ris_las$P, inputs), error = function(e) NULL)
      if (!is.null(q_las)) {
        dbh_las <- dbh_robusto(q_las)
        if (!is.na(dbh_las)) {
          risultato$DBHcyl_finale <- dbh_las
          risultato$metodo <- "ripopolamento_las_multifascia"
          risultato$n_fasce_valide <- ris_las$n_fasce_valide
          return(risultato)
        }
      }
    }
    
    P_taglio05 <- P_grezza[P_grezza[,3] >= 0.5, , drop = FALSE]
    q_t05 <- tryCatch(treeqsm(P_taglio05, inputs_vecchi), error = function(e) NULL)
    if (!is.null(q_t05)) {
      dbh_t05 <- dbh_robusto(q_t05)
      if (!is.na(dbh_t05)) { risultato$DBHcyl_finale <- dbh_t05; risultato$metodo <- "taglio05_vecchi_param" }
      else { risultato$metodo <- "taglio05_fallito_tutti_degeneri" }
    }
    due_fusti <- tryCatch(traccia_due_alberi_da_altezza(P_grezza), error = function(e) NULL)
    if (!is.null(due_fusti)) {
      A_t <- due_fusti$albero_A[due_fusti$albero_A[,3] >= 0.5, , drop = FALSE]
      B_t <- due_fusti$albero_B[due_fusti$albero_B[,3] >= 0.5, , drop = FALSE]
      qA <- tryCatch(treeqsm(A_t, inputs_vecchi), error = function(e) NULL)
      qB <- tryCatch(treeqsm(B_t, inputs_vecchi), error = function(e) NULL)
      if (!is.null(qA) && !is.null(qB)) {
        dA <- dbh_robusto(qA); dB <- dbh_robusto(qB)
        if (!is.na(dA) && !is.na(dB) && min(dA,dB)/max(dA,dB) > 0.75) {
          risultato$DBHcyl_finale <- sqrt(dA^2+dB^2); risultato$metodo <- "diametro_equivalente"
        }
      }
    }
    risultato$avviso <- "Diagnosi di contaminazione localizzata; il ripopolamento dal LAS non ha prodotto un risultato sufficientemente affidabile (segnale multifascia assente o incoerente). La stima riportata usa il metodo di fallback e potrebbe essere affetta da incertezza superiore alla norma."
    if (!is.null(ris_las)) risultato$motivo_las <- if (isTRUE(ris_las$successo)) paste("LAS letto ma non affidabile (n_fasce=", ris_las$n_fasce_valide, ", raggio=", round(ris_las$raggio_sicuro,2), ")") else ris_las$motivo
    return(risultato)
  }
  
  if (e_catena) {
    is_tronco0 <- cyl$branch == 1
    r_c0 <- cyl$radius[is_tronco0] * 100
    km_cat <- tryCatch(kmeans(r_c0, centers = 2), error = function(e) NULL)
    if (!is.null(km_cat)) {
      cluster_alto0 <- which.max(km_cat$centers)
      risultato$DBHcyl_finale <- mean(r_c0[km_cat$cluster == cluster_alto0]) * 2
      risultato$diagnosi <- "catena_discontinua"; risultato$metodo <- "clustering_raggio_catena"
      return(risultato)
    }
  }
  
  inputs_flat <- appiattisci_inputs(inputs)
  ref <- robust_trunk_reference(cyl, cyl$branch==1, seq_along(cyl$radius))
  Z_start <- ref$start[3]
  cover1 <- cover_sets(P_pulito, inputs_flat)
  ts1 <- tree_sets(P_pulito, cover1, inputs_flat)
  seg1 <- correct_segments(P_pulito, ts1$cover, segments(ts1$cover, ts1$Base, ts1$Forb), inputs_flat,0,1,1)
  RS <- relative_size(P_pulito, ts1$cover, seg1)
  cover2 <- cover_sets(P_pulito, inputs_flat, RS)
  ts2 <- tree_sets(P_pulito, cover2, inputs_flat, seg1)
  seg2 <- correct_segments(P_pulito, ts2$cover, segments(ts2$cover, ts2$Base, ts2$Forb), inputs_flat,1,1,0)
  trunk_idx <- unique(unlist(ts2$cover$ball[unlist(seg2$segments[[1]])]))
  trunk_pts <- P_pulito[trunk_idx, , drop = FALSE]
  V <- sweep(trunk_pts, 2, ref$start); hh <- as.vector(V %*% ref$axis)
  fasce_assolute <- list(c(1.2,1.4), c(1.1,1.6), c(0.9,1.4))
  if (Z_start > 1.35) fasce_assolute <- c(fasce_assolute, list(c(Z_start+0.02, Z_start+0.42)))
  migliore_sd <- Inf
  for (fascia_ass in fasce_assolute) {
    fascia_hh <- fascia_ass - Z_start
    if (fascia_hh[2] <= 0.05) next
    mask_zona <- hh > fascia_hh[1] & hh < fascia_hh[2]
    pts_zona <- trunk_pts[mask_zona, , drop = FALSE]
    if (nrow(pts_zona) < 30) next
    d_zona <- distances_to_line(pts_zona, ref$axis, ref$start)$d * 100
    mediana_grezza <- median(d_zona)
    if (mediana_grezza > raggio_max_plausibile) next
    idx_correnti <- seq_along(d_zona); rapporto_prec <- NA; idx_migliori <- idx_correnti
    for (iter in 1:max_iter) {
      km <- tryCatch(kmeans(d_zona[idx_correnti], centers = 2), error = function(e) NULL)
      if (is.null(km)) break
      centri <- sort(km$centers); rapporto <- centri[2]/centri[1]
      cluster_int <- which.min(km$centers); idx_nuovi <- idx_correnti[km$cluster == cluster_int]
      calo <- if(!is.na(rapporto_prec)) (rapporto_prec - rapporto)/rapporto_prec else NA
      if (!is.na(calo) && calo > soglia_calo) break
      if (length(idx_nuovi) < 20) break
      idx_migliori <- idx_nuovi; idx_correnti <- idx_nuovi; rapporto_prec <- rapporto
    }
    d_migliori <- d_zona[idx_migliori]
    if (mean(d_migliori) < soglia_restringimento * mediana_grezza) next
    pts_finali <- pts_zona[idx_migliori, , drop = FALSE]
    cyl0 <- list(radius = mean(d_migliori)/100, start = ref$start, axis = ref$axis)
    cyl_fit <- tryCatch(least_squares_cylinder(pts_finali, cyl0), error = function(e) NULL)
    if (is.null(cyl_fit) || !isTRUE(cyl_fit$conv)) next
    d_residui <- distances_to_line(pts_finali, cyl_fit$axis, cyl_fit$start)$d * 100
    sd_relativo <- sd(d_residui - cyl_fit$radius*100) / (cyl_fit$radius*100)
    if (sd_relativo < migliore_sd && sd_relativo <= soglia_residui) {
      migliore_sd <- sd_relativo
      risultato$DBHcyl_finale <- mean(d_migliori) * 2
      risultato$metodo <- "clustering_universale"
      risultato$sd_relativo <- round(sd_relativo,3)
      risultato$fascia_usata <- paste(round(fascia_ass,2), collapse="-")
    }
  }
  risultato
}

# --- PIPELINE FINALE v3 (Day 22, seconda parte): criterio migliorato per --
# clustering_universale. Scoperto (caso: albero 3) che scegliere "la fascia
# con sd_relativo minimo" tra quelle valide favorisce sistematicamente
# fit su campioni piccoli (piu' facili da adattare bene, sd basso) invece
# del fit piu' rappresentativo (piu' punti, sd leggermente piu' alto ma
# sotto soglia). Corretto: raccogliere i DBH di TUTTE le fasce che passano
# la soglia (ora 0.08, non piu' 0.075 - verificato che 0.075 lasciava
# spesso una sola fascia valida, rendendo il fix inefficace) e usarne la
# MEDIANA, non il minimo sd. Verificato: albero 3 migliora da 35.2% a 27.1%,
# albero 14 (stesso ramo) resta invariato a 1.1% - nessun danno collaterale.
# Risultato sui 15 alberi con risultato: mediana 8.8% (invariata), media
# 17.69% (da 18.23%).

pipeline_integrata_v11 <- function(P_grezza_input, inputs, inputs_vecchi, n_modelli = 5,
                                   soglia_calo = 0.07, soglia_residui = 0.08,
                                   max_iter = 6, raggio_max_plausibile = 80,
                                   soglia_spostamento_centro = 0.7, soglia_ricentraggio = 0.5,
                                   soglia_restringimento = 0.5, cx_nominale = NULL, cy_nominale = NULL) {
  
  P_grezza <- P_grezza_input
  disp <- rileva_dispersione_anomala(P_grezza)
  P_pulito <- P_grezza[P_grezza[,3] >= disp$h_inizio_tronco_pulito, , drop = FALSE]
  q <- treeqsm(P_pulito, inputs)
  cyl <- modello_non_degenere(q)$cylinder
  
  is_tronco_r <- cyl$branch == 1
  x_tronco <- median(cyl$start[is_tronco_r, 1]); y_tronco <- median(cyl$start[is_tronco_r, 2])
  spostamento_iniziale <- sqrt(x_tronco^2 + y_tronco^2)
  ricentrato <- FALSE
  if (spostamento_iniziale > soglia_ricentraggio) {
    P_grezza[,1] <- P_grezza[,1] - x_tronco; P_grezza[,2] <- P_grezza[,2] - y_tronco
    P_pulito[,1] <- P_pulito[,1] - x_tronco; P_pulito[,2] <- P_pulito[,2] - y_tronco
    ricentrato <- TRUE
    q <- treeqsm(P_pulito, inputs)
    cyl <- modello_non_degenere(q)$cylinder
  }
  
  dbh_standard <- dbh_robusto(q)
  if (is.na(dbh_standard)) dbh_standard <- q[[1]]$treedata$DBHcyl * 100
  
  diag <- diagnostica_rumore_interno_diretta(P_grezza)
  e_catena <- diagnostica_catena_discontinua(cyl)
  
  ref_base <- robust_trunk_reference(cyl, cyl$branch==1, seq_along(cyl$radius))
  dist_centro_iniziale <- sqrt(ref_base$start[1]^2 + ref_base$start[2]^2)
  usa_fascia_bassa <- dist_centro_iniziale > soglia_spostamento_centro
  
  applica_knn <- e_catena || (!is.null(diag$esito) && diag$esito == "RUMORE_INTERNO_DIFFUSO" && !usa_fascia_bassa)
  
  if (applica_knn) {
    keep_knn <- tryCatch(filtro_knn_fedele(P_grezza), error=function(e) rep(TRUE, nrow(P_grezza)))
    P_grezza_knn <- P_grezza[keep_knn, , drop = FALSE]
    disp2 <- tryCatch(rileva_dispersione_anomala(P_grezza_knn), error=function(e) NULL)
    if (!is.null(disp2)) {
      P_pulito2 <- P_grezza_knn[P_grezza_knn[,3] >= disp2$h_inizio_tronco_pulito, , drop = FALSE]
      q2 <- tryCatch(treeqsm(P_pulito2, inputs), error=function(e) NULL)
      if (!is.null(q2)) {
        P_grezza <- P_grezza_knn; P_pulito <- P_pulito2
        cyl <- modello_non_degenere(q2)$cylinder
        dbh_standard2 <- dbh_robusto(q2)
        dbh_standard <- if (is.na(dbh_standard2)) q2[[1]]$treedata$DBHcyl * 100 else dbh_standard2
      }
    }
  }
  
  risultato <- list(DBHcyl_finale = dbh_standard, DBHcyl_standard = dbh_standard,
                    diagnosi = diag$esito, metodo = "standard", ricentrato = ricentrato, avviso = NULL)
  
  if (!is.null(diag$esito) && diag$esito == "RUMORE_INTERNO_DIFFUSO") {
    dbh_corretti <- c()
    for (m in seq_len(n_modelli)) {
      inputs_m <- inputs; inputs_m$model <- m
      q_m <- tryCatch(treeqsm(P_pulito, inputs_m), error = function(e) NULL)
      if (is.null(q_m)) next
      inputs_m_flat <- appiattisci_inputs(inputs_m)
      cyl_m <- modello_non_degenere(q_m)$cylinder
      ref_m <- robust_trunk_reference(cyl_m, cyl_m$branch==1, seq_along(cyl_m$radius))
      Z_start_m <- ref_m$start[3]
      cover1_m <- cover_sets(P_pulito, inputs_m_flat)
      ts1_m <- tree_sets(P_pulito, cover1_m, inputs_m_flat)
      seg1_m <- correct_segments(P_pulito, ts1_m$cover, segments(ts1_m$cover, ts1_m$Base, ts1_m$Forb), inputs_m_flat,0,1,1)
      RS_m <- relative_size(P_pulito, ts1_m$cover, seg1_m)
      cover2_m <- cover_sets(P_pulito, inputs_m_flat, RS_m)
      ts2_m <- tree_sets(P_pulito, cover2_m, inputs_m_flat, seg1_m)
      seg2_m <- correct_segments(P_pulito, ts2_m$cover, segments(ts2_m$cover, ts2_m$Base, ts2_m$Forb), inputs_m_flat,1,1,0)
      trunk_idx_m <- unique(unlist(ts2_m$cover$ball[unlist(seg2_m$segments[[1]])]))
      trunk_pts_m <- P_pulito[trunk_idx_m, , drop = FALSE]
      V <- sweep(trunk_pts_m,2,ref_m$start); hh_m <- as.vector(V %*% ref_m$axis)
      if (Z_start_m > 1.35) { fascia_ass <- c(Z_start_m+0.02, Z_start_m+0.42) }
      else if (usa_fascia_bassa) { fascia_ass <- c(1.2, 1.4) } else { fascia_ass <- c(1.7, 2.1) }
      fascia_hh <- fascia_ass - Z_start_m
      if (fascia_hh[2] <= 0.05) next
      Tpts <- trunk_pts_m[hh_m > fascia_hh[1] & hh_m < fascia_hh[2], , drop = FALSE]
      if (nrow(Tpts) < 50) next
      dl_init <- distances_to_line(Tpts, ref_m$axis, ref_m$start)
      km <- tryCatch(kmeans(dl_init$d, centers = 2), error = function(e) NULL)
      if (is.null(km)) next
      mask_ext <- km$cluster == which.max(km$centers)
      dbh_corretti <- c(dbh_corretti, mean(dl_init$d[mask_ext])*2*100)
    }
    if (length(dbh_corretti) > 0) { risultato$DBHcyl_finale <- median(dbh_corretti); risultato$metodo <- "multimodello_centroguidato" }
    return(risultato)
  }
  
  if (!is.null(diag$esito) && diag$esito == "contaminazione_localizzata") {
    ris_las <- NULL
    if (!is.null(cx_nominale) && !is.null(cy_nominale)) {
      ris_las <- tryCatch(arricchisci_da_las(cx_nominale, cy_nominale), error = function(e) list(successo = FALSE, motivo = conditionMessage(e)))
    }
    if (!is.null(ris_las) && isTRUE(ris_las$successo) && isTRUE(ris_las$affidabile)) {
      q_las <- tryCatch(treeqsm(ris_las$P, inputs), error = function(e) NULL)
      if (!is.null(q_las)) {
        dbh_las <- dbh_robusto(q_las)
        if (!is.na(dbh_las)) {
          risultato$DBHcyl_finale <- dbh_las
          risultato$metodo <- "ripopolamento_las_multifascia"
          risultato$n_fasce_valide <- ris_las$n_fasce_valide
          return(risultato)
        }
      }
    }
    P_taglio05 <- P_grezza[P_grezza[,3] >= 0.5, , drop = FALSE]
    q_t05 <- tryCatch(treeqsm(P_taglio05, inputs_vecchi), error = function(e) NULL)
    if (!is.null(q_t05)) {
      dbh_t05 <- dbh_robusto(q_t05)
      if (!is.na(dbh_t05)) { risultato$DBHcyl_finale <- dbh_t05; risultato$metodo <- "taglio05_vecchi_param" }
      else { risultato$metodo <- "taglio05_fallito_tutti_degeneri" }
    }
    due_fusti <- tryCatch(traccia_due_alberi_da_altezza(P_grezza), error = function(e) NULL)
    if (!is.null(due_fusti)) {
      A_t <- due_fusti$albero_A[due_fusti$albero_A[,3] >= 0.5, , drop = FALSE]
      B_t <- due_fusti$albero_B[due_fusti$albero_B[,3] >= 0.5, , drop = FALSE]
      qA <- tryCatch(treeqsm(A_t, inputs_vecchi), error = function(e) NULL)
      qB <- tryCatch(treeqsm(B_t, inputs_vecchi), error = function(e) NULL)
      if (!is.null(qA) && !is.null(qB)) {
        dA <- dbh_robusto(qA); dB <- dbh_robusto(qB)
        if (!is.na(dA) && !is.na(dB) && min(dA,dB)/max(dA,dB) > 0.75) {
          risultato$DBHcyl_finale <- sqrt(dA^2+dB^2); risultato$metodo <- "diametro_equivalente"
        }
      }
    }
    risultato$avviso <- "Diagnosi di contaminazione localizzata; il ripopolamento dal LAS non ha prodotto un risultato sufficientemente affidabile."
    return(risultato)
  }
  
  if (e_catena) {
    is_tronco0 <- cyl$branch == 1
    r_c0 <- cyl$radius[is_tronco0] * 100
    km_cat <- tryCatch(kmeans(r_c0, centers = 2), error = function(e) NULL)
    if (!is.null(km_cat)) {
      cluster_alto0 <- which.max(km_cat$centers)
      risultato$DBHcyl_finale <- mean(r_c0[km_cat$cluster == cluster_alto0]) * 2
      risultato$diagnosi <- "catena_discontinua"; risultato$metodo <- "clustering_raggio_catena"
      return(risultato)
    }
  }
  
  inputs_flat <- appiattisci_inputs(inputs)
  ref <- robust_trunk_reference(cyl, cyl$branch==1, seq_along(cyl$radius))
  Z_start <- ref$start[3]
  cover1 <- cover_sets(P_pulito, inputs_flat)
  ts1 <- tree_sets(P_pulito, cover1, inputs_flat)
  seg1 <- correct_segments(P_pulito, ts1$cover, segments(ts1$cover, ts1$Base, ts1$Forb), inputs_flat,0,1,1)
  RS <- relative_size(P_pulito, ts1$cover, seg1)
  cover2 <- cover_sets(P_pulito, inputs_flat, RS)
  ts2 <- tree_sets(P_pulito, cover2, inputs_flat, seg1)
  seg2 <- correct_segments(P_pulito, ts2$cover, segments(ts2$cover, ts2$Base, ts2$Forb), inputs_flat,1,1,0)
  trunk_idx <- unique(unlist(ts2$cover$ball[unlist(seg2$segments[[1]])]))
  trunk_pts <- P_pulito[trunk_idx, , drop = FALSE]
  V <- sweep(trunk_pts, 2, ref$start); hh <- as.vector(V %*% ref$axis)
  fasce_assolute <- list(c(1.2,1.4), c(1.1,1.6), c(0.9,1.4))
  if (Z_start > 1.35) fasce_assolute <- c(fasce_assolute, list(c(Z_start+0.02, Z_start+0.42)))
  
  dbh_candidati <- c()
  for (fascia_ass in fasce_assolute) {
    fascia_hh <- fascia_ass - Z_start
    if (fascia_hh[2] <= 0.05) next
    mask_zona <- hh > fascia_hh[1] & hh < fascia_hh[2]
    pts_zona <- trunk_pts[mask_zona, , drop = FALSE]
    if (nrow(pts_zona) < 30) next
    d_zona <- distances_to_line(pts_zona, ref$axis, ref$start)$d * 100
    mediana_grezza <- median(d_zona)
    if (mediana_grezza > raggio_max_plausibile) next
    idx_correnti <- seq_along(d_zona); rapporto_prec <- NA; idx_migliori <- idx_correnti
    for (iter in 1:max_iter) {
      km <- tryCatch(kmeans(d_zona[idx_correnti], centers = 2), error = function(e) NULL)
      if (is.null(km)) break
      centri <- sort(km$centers); rapporto <- centri[2]/centri[1]
      cluster_int <- which.min(km$centers); idx_nuovi <- idx_correnti[km$cluster == cluster_int]
      calo <- if(!is.na(rapporto_prec)) (rapporto_prec - rapporto)/rapporto_prec else NA
      if (!is.na(calo) && calo > soglia_calo) break
      if (length(idx_nuovi) < 20) break
      idx_migliori <- idx_nuovi; idx_correnti <- idx_nuovi; rapporto_prec <- rapporto
    }
    d_migliori <- d_zona[idx_migliori]
    if (mean(d_migliori) < soglia_restringimento * mediana_grezza) next
    pts_finali <- pts_zona[idx_migliori, , drop = FALSE]
    cyl0 <- list(radius = mean(d_migliori)/100, start = ref$start, axis = ref$axis)
    cyl_fit <- tryCatch(least_squares_cylinder(pts_finali, cyl0), error = function(e) NULL)
    if (is.null(cyl_fit) || !isTRUE(cyl_fit$conv)) next
    d_residui <- distances_to_line(pts_finali, cyl_fit$axis, cyl_fit$start)$d * 100
    sd_relativo <- sd(d_residui - cyl_fit$radius*100) / (cyl_fit$radius*100)
    if (sd_relativo <= soglia_residui) {
      dbh_candidati <- c(dbh_candidati, mean(d_migliori) * 2)
    }
  }
  if (length(dbh_candidati) > 0) {
    risultato$DBHcyl_finale <- median(dbh_candidati)
    risultato$metodo <- "clustering_universale"
    risultato$n_fasce_usate <- length(dbh_candidati)
  }
  risultato
}
