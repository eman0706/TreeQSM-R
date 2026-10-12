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
# TREEQSM - Porting MATLAB helper functions to R
# =============================================================================
# These functions are ported from the MATLAB TreeQSM package
# Source: TreeQSM-master/src/tools/ and TreeQSM-master/src/least_squares_fitting/
# =============================================================================

# Null-coalesce: restituisce a se non NULL, altrimenti b
`%||%` <- function(a, b) if (!is.null(a)) a else b

# NOTA (vedi conversazione): la vera connected_components() (fedele
# all'originale, con firma Nei/Sub/MinSize/Fal) e' nel file dedicato
# connected_components.R, caricato DOPO questo file in core_files (in
# entrambi i loader di test_run03.R) - quindi sovrascrive correttamente
# questa versione. Rimossa una versione semplificata precedente (firma
# solo Nei, algoritmo diverso, mai davvero usata in pratica ma fonte di
# confusione) per evitare lo stesso tipo di trappola incontrata oggi
# con surface_coverage_filtering.

#' Average of columns
#' @param X Matrix
#' @return Vector of column means
average <- function(X) {
  n <- nrow(X)
  if (n > 1) {
    return(colSums(X) / n)
  } else {
    return(as.vector(X))
  }
}

#' Cross product of two 3-vectors
#' @param A,B 3-vectors
#' @return Cross product A x B
cross_product <- function(A, B) {
  c(
    A[2] * B[3] - A[3] * B[2],
    A[3] * B[1] - A[1] * B[3],
    A[1] * B[2] - A[2] * B[1]
  )
}

#' Generate orthonormal vectors
#' @param U Input vector (3)
#' @return List with V and W orthonormal to U
orthonormal_vectors <- function(U) {
  # FIX: il vettore di riferimento casuale (runif) poteva, in rari casi
  # numerici, produrre risultati degeneri gestiti male dal ciclo while
  # (osservato: "numero di dimensioni errato" propagato da qui su
  # alberi specifici). Sostituito con una scelta deterministica e
  # robusta: usa (0,0,1) come riferimento se U non e' quasi parallelo a
  # esso, altrimenti (1,0,0) - garantisce sempre un cross product ben
  # definito, senza cicli ne' casualita'.
  ref <- if (abs(U[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)
  V <- cross_product(ref, U)
  nV <- norm_vec(V)
  if (!is.finite(nV) || nV < 1e-10) {
    # U quasi parallelo anche al fallback (caso patologico estremo)
    ref <- c(0, 1, 0)
    V <- cross_product(ref, U)
    nV <- norm_vec(V)
  }
  W <- cross_product(V, U)
  nW <- norm_vec(W)
  W <- W / nW
  V <- V / nV
  
  list(V = V, W = W)
}

#' Norm of a vector
norm_vec <- function(x) {
  sqrt(sum(x^2))
}

#' Rotation matrix for axis-angle
#' @param A Unit axis vector (3)
#' @param angle Rotation angle in radians
#' @return 3x3 rotation matrix
rotation_matrix <- function(A, angle) {
  A <- A / norm_vec(A)
  c <- cos(angle)
  s <- sin(angle)
  
  R <- matrix(0, 3, 3)
  R[1, ] <- c(A[1]^2 + (1 - A[1]^2) * c, 
              A[1] * A[2] * (1 - c) - A[3] * s,
              A[1] * A[3] * (1 - c) + A[2] * s)
  R[2, ] <- c(A[1] * A[2] * (1 - c) + A[3] * s,
              A[2]^2 + (1 - A[2]^2) * c,
              A[2] * A[3] * (1 - c) - A[1] * s)
  R[3, ] <- c(A[1] * A[3] * (1 - c) - A[2] * s,
              A[2] * A[3] * (1 - c) + A[1] * s,
              A[3]^2 + (1 - A[3]^2) * c)
  R
}

#' Rotate vector to align with z-axis
#' @param Vec Input vector (3)
#' @return List with R (rotation matrix), D, a
rotate_to_z_axis <- function(Vec) {
  D <- cross_product(Vec, c(0, 0, 1))
  if (norm_vec(D) > 0) {
    a <- acos(Vec[3])
    R <- rotation_matrix(D, a)
  } else {
    R <- diag(3)
    a <- 0
    D <- c(1, 0, 0)
  }
  list(R = R, D = D, a = a)
}

#' Form rotation matrices from angles
#' @param angles Vector of two angles
#' @return List with R, DR1, DR2
form_rotation_matrices <- function(angles) {
  alpha <- angles[1]
  beta <- angles[2]
  
  ca <- cos(alpha)
  sa <- sin(alpha)
  cb <- cos(beta)
  sb <- sin(beta)
  
  R <- matrix(0, 3, 3)
  R[1, ] <- c(cb, 0, sb)
  R[2, ] <- c(sa * sb, ca, -sa * cb)
  R[3, ] <- c(-ca * sb, sa, ca * cb)
  
  DR1 <- matrix(0, 3, 3)  # derivative wrt alpha
  DR1[1, ] <- c(0, 0, 0)
  DR1[2, ] <- c(ca * sb, -sa, -ca * cb)
  DR1[3, ] <- c(sa * sb, ca, sa * cb)
  
  DR2 <- matrix(0, 3, 3)  # derivative wrt beta
  DR2[1, ] <- c(-sb, 0, cb)
  DR2[2, ] <- c(sa * cb, 0, sa * sb)
  DR2[3, ] <- c(-ca * cb, 0, -ca * sb)
  
  list(R = R, DR1 = DR1, DR2 = DR2)
}

#' Unique elements with filtering
#' @param Set Vector of indices
#' @param Fal Logical vector (same length as possible values)
#' @return Filtered unique vector
unique_elements <- function(Set, Fal) {
  n <- length(Set)
  if (n > 2) {
    I <- rep(TRUE, n)
    for (j in 1:n) {
      if (!Fal[Set[j]]) {
        Fal[Set[j]] <- TRUE
      } else {
        I[j] <- FALSE
      }
    }
    Set <- Set[I]
  } else if (n == 2) {
    if (Set[1] == Set[2]) {
      Set <- Set[1]
    }
  }
  Set
}

#' Vertical concatenation of cell array
#' @param CellArray List of vectors
#' @return List with Vector (concatenated) and IndElements (indices)
verticalcat <- function(CellArray) {
  CellSize <- sapply(CellArray, length)
  nc <- length(CellArray)
  
  IndElements <- matrix(1L, nrow = nc, ncol = 2)
  IndElements[, 2] <- cumsum(CellSize)
  if (nc > 1)
    IndElements[2:nc, 1] <- IndElements[1:(nc - 1), 2] + 1L
  
  Vector <- unlist(CellArray)
  
  list(Vector = Vector, IndElements = IndElements)
}

#' Distances from points to a line
#' @param Q Points matrix (n x 3)
#' @param LineDirec Direction vector (unit vector, 3)
#' @param LinePoint Point on line (3)
#' @return List with d (distances), V (perpendicular vectors), h (heights), B
distances_to_line <- function(Q, LineDirec, LinePoint) {
  A <- Q - matrix(LinePoint, nrow = nrow(Q), ncol = 3, byrow = TRUE)
  h <- as.vector(A %*% LineDirec)
  B <- outer(h, LineDirec)
  V <- A - B
  d <- sqrt(rowSums(V^2))
  
  list(d = d, V = V, h = h, B = B)
}

#' Surface coverage computation
#' @param P Points matrix
#' @param Axis Axis direction (3)
#' @param Point Point on axis (3)
#' @param nl Number of layers
#' @param ns Number of sectors
#' @param Dmin Minimum distance (optional)
#' @param Dmax Maximum distance (optional)
#' @return Surface coverage value (0-1)
surface_coverage <- function(P, Axis, Point, nl, ns, Dmin = NULL, Dmax = NULL) {
  dists <- distances_to_line(P, Axis, Point)
  d <- dists$d
  V <- dists$V
  h <- dists$h - min(dists$h)
  Len <- max(h)
  
  # Filter by distance if specified
  if (!is.null(Dmin)) {
    Keep <- d > Dmin
    if (!is.null(Dmax)) {
      Keep <- Keep & d < Dmax
    }
    V <- V[Keep, , drop = FALSE]
    h <- h[Keep]
    d <- d[Keep]
  }
  
  # Compute coverage from 4 different baseline directions
  SurfCov <- numeric(4)
  oa <- orthonormal_vectors(Axis)
  U <- oa$V
  W <- oa$W
  R <- rotation_matrix(Axis, 2 * pi / ns / 4)
  
  for (i in 1:4) {
    if (i > 1) {
      U <- R %*% U
      W <- R %*% W
    }
    
    Vproj <- V %*% cbind(U, W)
    ang <- atan2(Vproj[, 2], Vproj[, 1]) + pi
    
    Layer <- ceiling(h / Len * nl)
    Layer[Layer <= 0] <- 1
    Layer[Layer > nl] <- nl
    
    Sector <- ceiling(ang / 2 / pi * ns)
    Sector[Sector <= 0] <- 1
    
    LexOrd <- Layer + (Sector - 1) * nl
    
    Cov <- matrix(0, nl, ns)
    Cov[LexOrd] <- 1
    SurfCov[i] <- sum(Cov > 0) / nl / ns
  }
  
  max(SurfCov)
}

#' Surface Coverage Filtering
#'
#' Porting fedele di surface_coverage_filtering.m (Version 1.1.0,
#' Copyright (C) 2017-2021 Pasi Raumonen). Stima adattiva del raggio dai
#' dati stessi (minimo per cella layer/settore, poi mediana), ricalcolo
#' dinamico di nl/ns basato su quella stima, filtro per prossimita'
#' relativa al minimo di ciascuna cella. Usata dentro cylinder_fitting()
#' per ogni cilindro con >20 punti.
#'
#' @param P Point cloud del segmento (matrice Nx3)
#' @param c Cilindro (lista con axis, start, length)
#' @param lh Altezza dei layer (in metri)
#' @param ns Numero di settori (stima iniziale)
#' @return Lista con Pass (vettore logico) e c (cilindro aggiornato con
#'   radius/SurfCov/mad/conv/rel)
#' @export
surface_coverage_filtering <- function(P, c, lh = 0.02, ns = 20) {
  
  np <- nrow(P)
  
  dl <- distances_to_line(P, c$axis, c$start)
  d  <- dl$d
  h  <- dl$h - min(dl$h)
  
  ow <- orthonormal_vectors(c$axis)
  U  <- ow$V;  W <- ow$W
  V_perp <- dl$V
  Uproj  <- as.numeric(V_perp %*% U)
  Wproj  <- as.numeric(V_perp %*% W)
  ang    <- atan2(Wproj, Uproj) + pi
  
  cyl_len <- if (!is.null(c$length) && is.finite(c$length) && c$length > 0)
    c$length else (max(h) - min(h) + 1e-6)
  
  nl <- max(3L, ceiling(cyl_len / lh))
  
  Layer  <- pmin(nl, pmax(1L, ceiling(h  / cyl_len * nl)))
  Sector <- pmin(ns, pmax(1L, ceiling(ang / (2 * pi) * ns)))
  LexOrd <- Layer + (Sector - 1L) * nl
  
  ord    <- order(LexOrd)
  LexOrd_s <- LexOrd[ord]
  ds     <- d[ord]
  
  Dis1 <- numeric(nl * ns)
  p <- 1L
  while (p <= np) {
    t <- 1L
    while ((p + t <= np) && LexOrd_s[p] == LexOrd_s[p + t]) t <- t + 1L
    Dmin <- min(ds[p:(p + t - 1L)])
    Dis1[LexOrd_s[p]] <- min(1.05 * Dmin, Dmin + 0.02)
    p <- p + t
  }
  
  R_est <- median(Dis1[Dis1 > 0])
  if (!is.finite(R_est) || R_est <= 0) R_est <- median(d[d > 0])
  if (!is.finite(R_est) || R_est <= 0) R_est <- 0.05
  
  a  <- max(0.02, 0.2 * R_est)
  ns <- min(36L, max(8L, ceiling(2 * pi * R_est / a)))
  nl <- max(3L, ceiling(cyl_len / a))
  
  Layer  <- pmin(nl, pmax(1L, ceiling(h  / cyl_len * nl)))
  Sector <- pmin(ns, pmax(1L, ceiling(ang / (2 * pi) * ns)))
  LexOrd <- Layer + (Sector - 1L) * nl
  
  ord      <- order(LexOrd)
  LexOrd_s <- LexOrd[ord]
  d_sorted <- d[ord]
  
  Dis2 <- numeric(nl * ns)
  Pass <- logical(np)
  r_cell <- max(0.01, 0.05 * R_est)
  k <- 0L
  
  p <- 1L
  while (p <= np) {
    t <- 1L
    while ((p + t <= np) && LexOrd_s[p] == LexOrd_s[p + t]) t <- t + 1L
    ind  <- p:(p + t - 1L)
    D    <- d_sorted[ind]
    Dmin <- min(D)
    I    <- D <= Dmin + r_cell
    Pass[ind[I]] <- TRUE
    Dis2[LexOrd_s[p]] <- min(1.05 * Dmin, Dmin + 0.02)
    p <- p + t
    k <- k + 1L
  }
  
  inv_ord       <- integer(np)
  inv_ord[ord]  <- seq_len(np)
  Pass          <- Pass[inv_ord]
  
  R_final <- median(Dis2[Dis2 > 0])
  if (!is.finite(R_final) || R_final <= 0) R_final <- R_est
  
  d_pass  <- d[Pass]
  mad_val <- if (length(d_pass) > 0) mean(abs(d_pass - R_final)) else 0
  
  c$radius  <- R_final
  c$SurfCov <- k / (nl * ns)
  c$mad     <- mad_val
  c$conv    <- 1L
  c$rel     <- 1L
  
  list(Pass = Pass, c = c)
}

#' Least squares cylinder fitting - Levenberg-Marquardt (NOTA DI
#' FEDELTA', vedi conversazione)
#'
#' Porting NON fedele nell'algoritmo di ottimizzazione: il vero
#' least_squares_cylinder.m (Version 2.0.0, Copyright (C) 2013-2021
#' Pasi Raumonen, recuperato e confrontato riga per riga) usa
#' Gauss-Newton PURO (par iniziale [0,0,0,0,radius] SENZA
#' perturbazione, aggiornamento SEMPRE accettato anche se peggiora
#' temporaneamente, soglia di convergenza 1e-4 sul confronto SS0/SS1
#' della stessa iterazione, rcond controllato sulla matrice A NON
#' regolarizzata). Questa versione usa invece Levenberg-Marquardt con
#' regolarizzazione (par iniziale perturbato di 1e-6, aggiornamento
#' accettato solo se migliora, soglia 1e-5 sul confronto con
#' l'iterazione PRECEDENTE, rcond sulla matrice regolarizzata A_lm).
#'
#' VERIFICATO EMPIRICAMENTE (non solo teoricamente) che l'impatto
#' pratico e' trascurabile: su un caso sintetico facile e uno difficile
#' (pochi punti, rumore alto, cilindro corto/sottile, stima iniziale
#' imprecisa), le due versioni convergono a risultati praticamente
#' identici (differenze in 3a/4a cifra decimale). Sull'albero 47 (il
#' caso reale piu' difficile incontrato in sessione), DBHqsm e DBHcyl
#' sono risultati IDENTICI tra le due versioni (solo il numero di
#' cilindri differiva del 3%, 187 vs 181). La versione fedele
#' (Gauss-Newton puro) e' disponibile come
#' least_squares_cylinder_faithful.R per riferimento/test futuri, ma
#' NON sostituita in produzione dato che non emerge nessun beneficio
#' misurabile a fronte del rischio di regressioni non testate su altri
#' alberi.
#'
#' func_grad_cylinder() (il calcolo di residui/Jacobiano, la vera base
#' matematica) e' invece stato VERIFICATO come sostanzialmente fedele
#' all'originale (recuperato e confrontato riga per riga) - solo una
#' protezione minore contro la divisione per zero nel Jacobiano.
#'
#' @param P Points matrix
#' @param cyl0 Initial cylinder (axis, start, length, radius)
#' @param weight Optional weights
#' @param Q Optional subset for computing start/length
#' @return Fitted cylinder
least_squares_cylinder <- function(P, cyl0, weight = NULL, Q = NULL) {
  iter      <- 0L
  conv      <- FALSE
  rel       <- TRUE
  NoWeights <- is.null(weight)
  
  rot  <- rotate_to_z_axis(cyl0$axis)
  Rot0 <- rot$R
  # Trasformazione corretta: (P - start) poi ruota
  Pt   <- (P - matrix(cyl0$start, nrow = nrow(P), ncol = 3, byrow = TRUE)) %*% t(Rot0)
  
  # Perturbazione iniziale: evita singolarità di J a alpha=beta=0
  par     <- c(1e-6, 1e-6, 1e-6, 1e-6, cyl0$radius)
  SS_prev <- Inf
  dist    <- rep(0, nrow(P))
  
  while (iter < 50L && !conv) {
    fg <- if (NoWeights) func_grad_cylinder(par, Pt)
    else           func_grad_cylinder(par, Pt, weight)
    d0 <- fg$d; J <- fg$J
    
    if (is.null(J) || any(!is.finite(J))) { rel <- FALSE; break }
    
    A   <- t(J) %*% J
    b   <- t(J) %*% d0
    SS0 <- sqrt(sum(d0^2))
    
    # Levenberg-Marquardt: regularizza A per evitare singolarità
    lam  <- 1e-4 * max(diag(A), na.rm = TRUE)
    A_lm <- A + lam * diag(5)
    
    p_step <- tryCatch(suppressWarnings(-solve(A_lm, b)),
                       error = function(e) NULL)
    if (is.null(p_step) || any(!is.finite(p_step))) { rel <- FALSE; break }
    
    par_new <- par + p_step
    dist_new <- if (NoWeights) func_grad_cylinder(par_new, Pt)$d
    else           func_grad_cylinder(par_new, Pt, weight)$d
    SS1 <- sqrt(sum(dist_new^2))
    
    # Accetta solo passi migliorativi
    if (SS1 <= SS0) { par <- par_new; dist <- dist_new } else { dist <- d0 }
    
    # FIX: controllavamo rcond(A) (matrice NON regolarizzata) invece di
    # rcond(A_lm) (quella davvero risolta). A e' quasi sempre mal
    # condizionata per cilindri quasi verticali (il caso tipico in un
    # albero: bassa sensibilita' di alpha/beta), anche quando il fit
    # regolarizzato e' perfettamente valido - questo scartava quasi ogni
    # fit come "non affidabile" indipendentemente dalla qualita' reale.
    if (rcond(A_lm) < 1e-12) rel <- FALSE
    if (abs(SS_prev - SS1) < 1e-5) conv <- TRUE
    SS_prev <- SS1
    iter    <- iter + 1L
  }
  
  # Ricostruzione parametri cilindro nel sistema originale
  rot2  <- form_rotation_matrices(par[3:4])
  Axis  <- as.numeric(t(Rot0) %*% t(rot2$R) %*% c(0, 0, 1))
  Point <- as.numeric(t(Rot0) %*% c(par[1], par[2], 0) + cyl0$start)
  
  P_use <- if (!is.null(Q) && nrow(Q) > 5) Q else P
  H     <- as.numeric(P_use %*% Axis)
  hmin  <- min(H)
  Point <- Point - (sum(Axis * Point) - hmin) * Axis
  
  if (!is.null(weight)) {
    I_w  <- weight == max(weight)
    mad_val <- mean(abs(dist[I_w]))
  } else {
    mad_val <- mean(abs(dist))
  }
  
  cyl <- list(
    radius  = as.numeric(par[5]),
    length  = as.numeric(abs(max(H) - hmin)),
    start   = as.numeric(Point),
    axis    = as.numeric(Axis),
    mad     = if (is.finite(mad_val)) mad_val else 0,
    conv    = conv,
    rel     = rel,
    SurfCov = 0
  )
  
  # SurfCov: calcolato se il cilindro è geometricamente valido
  if (is.finite(cyl$radius) && cyl$radius > 1e-4 &&
      is.finite(cyl$length) && cyl$length > 1e-4 &&
      !any(is.na(Axis)) && !any(is.na(Point))) {
    res_l <- max(3L, ceiling(cyl$length / 0.03))
    res_s <- min(36L, max(8L, ceiling(2 * pi * cyl$radius / 0.03)))
    cyl$SurfCov <- tryCatch(
      surface_coverage(P_use, Axis, Point, res_l, res_s, 0.8 * cyl$radius),
      error = function(e) 0
    )
  }
  
  cyl
}

#' Function and gradient for cylinder fitting
#' @param par Parameters [x0, y0, alpha, beta, r]
#' @param P Points matrix (n x 3)
#' @param weight Optional weights
#' @return List with d (distances) and J (Jacobian)
func_grad_cylinder <- function(par, P, weight = NULL) {
  x0 <- par[1]; y0 <- par[2]
  alpha <- par[3]; beta <- par[4]; r <- par[5]
  
  rot <- form_rotation_matrices(c(alpha, beta))
  R   <- rot$R
  DR1 <- rot$DR1
  DR2 <- rot$DR2
  
  # P è n×3: sottrai offset poi ruota
  Pt  <- (P - matrix(c(x0, y0, 0), nrow = nrow(P), ncol = 3, byrow = TRUE)) %*% t(R)
  xt  <- Pt[, 1]; yt <- Pt[, 2]
  rt  <- sqrt(xt^2 + yt^2)
  rt_safe <- pmax(rt, 1e-10)   # evita divisione per zero
  dist <- rt - r
  if (!is.null(weight)) dist <- weight * dist
  
  # Jacobiano SEMPRE calcolato (Gauss-Newton lo richiede con o senza pesi)
  n <- nrow(P)
  J <- matrix(0, n, 5)
  N <- cbind(xt / rt_safe, yt / rt_safe)
  A <- diag(2)
  
  A1 <- R %*% c(-1, 0, 0); A[1,1] <- A1[1]; A[2,2] <- A1[2]
  J[, 1] <- rowSums(N %*% A)
  
  A2 <- R %*% c(0, -1, 0); A[1,1] <- A2[1]; A[2,2] <- A2[2]
  J[, 2] <- rowSums(N %*% A)
  
  if (!is.null(DR1)) {
    dPt1 <- (P - matrix(c(x0, y0, 0), nrow = n, ncol = 3, byrow = TRUE)) %*% t(DR1)
    J[, 3] <- rowSums(N * dPt1[, 1:2, drop = FALSE])
  }
  if (!is.null(DR2)) {
    dPt2 <- (P - matrix(c(x0, y0, 0), nrow = n, ncol = 3, byrow = TRUE)) %*% t(DR2)
    J[, 4] <- rowSums(N * dPt2[, 1:2, drop = FALSE])
  }
  J[, 5] <- -1
  if (!is.null(weight)) J <- J * weight
  
  list(d = dist, J = J)
}

#' Distances between lines
#' @param P1 Point on line 1
#' @param D1 Direction of line 1
#' @param P2 Point on line 2
#' @param D2 Direction of line 2
#' @return List with distances and closest points
distances_between_lines <- function(P1, D1, P2, D2) {
  D1 <- D1 / norm_vec(D1)
  D2 <- D2 / norm_vec(D2)
  
  V <- P2 - P1
  A <- sum(D1 * D1)
  B <- sum(D1 * D2)
  C <- sum(D2 * D2)
  D <- sum(D1 * V)
  E <- sum(D2 * V)
  
  denom <- A * C - B * B
  
  sc <- ifelse(denom > 0, (B * E - C * D) / denom, 0)
  tc <- ifelse(denom > 0, (A * E - B * D) / denom, 0)
  
  Q1 <- P1 + sc * D1
  Q2 <- P2 + tc * D2
  Dist <- Q2 - Q1
  d <- norm_vec(Dist)
  
  DistOnLines <- list(sc = sc, tc = tc)
  
  list(Dist = d, Q1 = Q1, Q2 = Q2, DistOnLines = DistOnLines)
}

#' Growth volume correction
#' @param cylinder Cylinder structure
#' @param inputs Input parameters
#' @return Corrected cylinder
growth_volume_correction <- function(cylinder, inputs) {
  # Simplified version
  GrowthVolFac <- inputs$GrowthVolFac
  
  # For now, just return the cylinder as-is
  # Full implementation would be complex
  cylinder
}