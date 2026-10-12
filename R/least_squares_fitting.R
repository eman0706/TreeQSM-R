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

# ==============================================================================
# LEAST_SQUARES_FITTING.R
# Porting ottimizzato per R con protezioni per singolarità e raggi nulli
# ==============================================================================

#' Calcola il prodotto vettoriale tra due vettori
cross_product <- function(v1, v2) {
  return(c(v1[2] * v2[3] - v1[3] * v2[2],
           v1[3] * v2[1] - v1[1] * v2[3],
           v1[1] * v2[2] - v1[2] * v2[1]))
}

#' Calcola la norma di un vettore
norm_vec <- function(v) {
  return(sqrt(sum(v^2)))
}

#' Genera una matrice di rotazione dato un asse e un angolo
rotation_matrix <- function(Axis, Angle) {
  Axis <- Axis / norm_vec(Axis)
  s <- sin(Angle)
  c <- cos(Angle)
  C <- 1 - c
  R <- matrix(0, 3, 3)
  R[1, 1] <- Axis[1]^2 * C + c
  R[1, 2] <- Axis[1] * Axis[2] * C - Axis[3] * s
  R[1, 3] <- Axis[1] * Axis[3] * C + Axis[2] * s
  R[2, 1] <- Axis[2] * Axis[1] * C + Axis[3] * s
  R[2, 2] <- Axis[2]^2 * C + c
  R[2, 3] <- Axis[2] * Axis[3] * C - Axis[1] * s
  R[3, 1] <- Axis[3] * Axis[1] * C - Axis[2] * s
  R[3, 2] <- Axis[3] * Axis[2] * C + Axis[1] * s
  R[3, 3] <- Axis[3]^2 * C + c
  return(R)
}

#' Ruota un vettore verso l'asse Z
rotate_to_z_axis <- function(Vec) {
  Vec <- Vec / norm_vec(Vec)
  # Se è già quasi allineato a Z
  if (abs(Vec[3]) > 0.9999) {
    return(list(R = diag(3), D = c(1, 0, 0), a = 0))
  }
  D <- cross_product(Vec, c(0, 0, 1))
  D <- D / norm_vec(D)
  a <- acos(Vec[3])
  R <- rotation_matrix(D, a)
  return(list(R = R, D = D, a = a))
}

#' Genera matrici di rotazione per piccole perturbazioni (alpha, beta)
form_rotation_matrices <- function(par) {
  a <- par[1]
  b <- par[2]
  # Matrice di rotazione approssimata per piccoli angoli
  R <- matrix(c(cos(b), 0, -sin(b),
                sin(a)*sin(b), cos(a), sin(a)*cos(b),
                cos(a)*sin(b), -sin(a), cos(a)*cos(b)), 3, 3, byrow = TRUE)
  return(list(R = R))
}

#' Funzione principale di fitting del cilindro
least_squares_cylinder <- function(P, cyl0, weight = NULL, Q = NULL) {
  # P: punti del segmento, cyl0: cilindro iniziale stimato
  
  maxiter <- 25
  iter <- 0
  conv <- FALSE
  rel <- TRUE
  
  # --- PROTEZIONE 1: Numero minimo di punti ---
  if (is.null(P) || nrow(P) < 3) {
    return(list(radius=0.005, length=0.1, start=cyl0$start, axis=cyl0$axis, conv=FALSE, rel=FALSE))
  }
  
  # Trasformazione: Allineamento all'asse Z iniziale
  rot_result <- rotate_to_z_axis(cyl0$axis)
  Rot0 <- rot_result$R
  Point0 <- cyl0$start
  Rad0 <- max(0.001, cyl0$radius) # Protezione raggio iniziale
  
  # Centra e ruota i punti
  Pt <- sweep(P, 2, Point0, "-") %*% t(Rot0)
  
  # Parametri iniziali per l'ottimizzazione (alpha, beta piccoli spostamenti)
  par <- c(0, 0) 
  
  # Gauss-Newton
  while (iter < maxiter && !conv) {
    # Calcola Jacobiano e distanze (usando la funzione ausiliaria func_grad_axis)
    # Nota: Assumiamo che func_grad_axis sia definita nel tuo ambiente
    fg <- func_grad_cylinder(par, Pt, Rad0) 
    d0 <- fg$dist
    J <- fg$J
    
    SS0 <- sum(d0^2)
    A <- t(J) %*% J
    b <- t(J) %*% d0
    
    # --- PROTEZIONE 2: Controllo singolarità matrice A ---
    if (rcond(A) < 1e-15) {
      rel <- FALSE
      break
    }
    
    # Risoluzione sistema lineare con tryCatch
    p <- tryCatch(
      solve(A, -b),
      error = function(e) NULL
    )
    
    if (is.null(p) || any(!is.finite(p))) {
      rel <- FALSE
      break
    }
    
    # Aggiornamento parametri
    par_new <- par + p
    
    # Calcolo nuovo errore
    fg_new <- func_grad_cylinder(par_new, Pt, Rad0)
    SS1 <- sum(fg_new$dist^2)
    
    # --- PROTEZIONE 3: Step-back se l'errore aumenta ---
    if (SS1 > SS0) {
      p <- 0.5 * p # Riduce il passo
      par <- par + p
    } else {
      par <- par_new
    }
    
    # Verifica convergenza
    if (abs(SS0 - SS1) / (SS0 + 1e-6) < 1e-5) {
      conv <- TRUE
    }
    
    iter <- iter + 1
  }
  
  # --- RICOSTRUZIONE GEOMETRIA FINALE ---
  rot_final <- form_rotation_matrices(par)
  # Nuovo asse ruotato nel sistema originale
  Axis_new <- as.vector(t(Rot0) %*% t(rot_final$R) %*% c(0, 0, 1))
  
  # Proiezione punti sull'asse per trovare lunghezza e start
  P_proj <- P %*% Axis_new
  h_min <- min(P_proj)
  h_max <- max(P_proj)
  L_final <- h_max - h_min
  
  # Nuovo punto di inizio (proiettato sulla base)
  # h0 è la posizione del vecchio start lungo il nuovo asse
  h0 <- sum(Point0 * Axis_new)
  Point_new <- Point0 + (h_min - h0) * Axis_new
  
  # --- PROTEZIONE 4: Fix raggio finale (il raggio è Rad0 + variazione se inclusa in par) ---
  # In questa implementazione semplificata usiamo il raggio calcolato dalla dispersione
  dist_final <- sqrt(rowSums((sweep(P, 2, Point_new, "-") %*% t(rotation_matrix(Axis_new, 0)))[,1:2]^2))
  Radius_final <- mean(dist_final)
  
  # Se il raggio calcolato è assurdo (0 o troppo grande), usa Rad0
  if (is.na(Radius_final) || Radius_final < 0.001) Radius_final <- Rad0
  if (Radius_final > 1.0) Radius_final <- Rad0 # Protezione alberi giganti/errori
  
  return(list(
    radius = as.numeric(Radius_final),
    length = as.numeric(L_final),
    start  = as.numeric(Point_new),
    axis   = as.numeric(Axis_new),
    conv   = conv,
    rel    = rel
  ))
}

#' Funzione di supporto per il gradiente (distanze dal cilindro)
func_grad_cylinder <- function(par, Pt, Rad) {
  # par = [alpha, beta]
  a <- par[1]
  b <- par[2]
  
  # Matrice di rotazione per i piccoli angoli
  # Approssimazione lineare delle rotazioni attorno a X e Y
  # x' = x*cos(b) + z*sin(b)
  # y' = x*sin(a)sin(b) + y*cos(a) - z*sin(a)cos(b)
  
  cos_a <- cos(a); sin_a <- sin(a)
  cos_b <- cos(b); sin_b <- sin(b)
  
  # Coordinate ruotate
  X <- Pt[,1]*cos_b + Pt[,3]*sin_b
  Y <- Pt[,1]*sin_a*sin_b + Pt[,2]*cos_a - Pt[,3]*sin_a*cos_b
  
  dist_sq <- X^2 + Y^2
  dist <- sqrt(dist_sq) - Rad
  
  # Jacobiano (semplificato per alpha e beta)
  # d_dist/d_alpha e d_dist/d_beta
  J <- matrix(0, nrow(Pt), 2)
  J[,1] <- Y / (dist + Rad) * (Pt[,1]*cos_a*sin_b - Pt[,2]*sin_a + Pt[,3]*cos_a*cos_b)
  J[,2] <- (X / (dist + Rad) * (-Pt[,1]*sin_b + Pt[,3]*cos_b)) + 
    (Y / (dist + Rad) * (Pt[,1]*sin_a*cos_b + Pt[,3]*sin_a*sin_b))
  
  return(list(dist = dist, J = J))
}