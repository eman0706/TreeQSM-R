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

#' Function and Gradient for Circle Centre
#'
#' NON presente nel sorgente recuperato (chiamata da
#' least_squares_circle_centre.m ma definita altrove nel repository
#' originale, mai recuperata) - implementazione propria, analoga a
#' func_grad_cylinder (gia' presente in helpers.R) ma per un cerchio 2D
#' a raggio FISSO (si stima solo il centro, non il raggio): residui
#' come distanza dal centro meno il raggio fisso, Jacobiano rispetto
#' alle sole coordinate del centro.
#'
#' @param P Punti 2D (matrice Nx2)
#' @param par Parametri correnti c(cx, cy, R) - R resta fisso durante
#'   l'iterazione, presente solo per uniformita' con l'interfaccia
#' @return Lista con dist (residui, vettore N) e J (Jacobiano, Nx2)
#' @export
func_grad_circle_centre <- function(P, par) {
  cx <- par[1]; cy <- par[2]; R <- par[3]
  dx <- P[, 1] - cx
  dy <- P[, 2] - cy
  r  <- sqrt(dx^2 + dy^2)
  dist <- r - R
  J <- cbind(-dx / r, -dy / r)
  list(dist = dist, J = J)
}

#' Least Squares Circle Centre
#'
#' Porting fedele di least_squares_circle_centre.m (Version 1.0.0,
#' Copyright (C) 2017-2021 Pasi Raumonen). Fit ai minimi quadrati di un
#' cerchio 2D con raggio DATO (fissato) - stima solo il centro, via
#' iterazione di Gauss-Newton.
#'
#' NOTA DI FEDELTA': l'originale, nel calcolo di ArcCov, non ha un
#' clamp esplicito per l'indice angolare quando risulta 0 dopo
#' ceil(ang/2/pi*100) (a differenza di altri punti del codice, es.
#' surface_coverage2.m, che clampano esplicitamente) - in MATLAB questo
#' indicizzerebbe con 0, dando errore. In R, I[0] <- TRUE non da
#' errore ma ignora silenziosamente quell'elemento (comportamento
#' diverso, ma piu' sicuro per la nostra pipeline - non introduce un
#' crash raro per un caso limite che l'originale stesso non gestisce
#' esplicitamente).
#'
#' @param P Punti 2D (matrice Nx2)
#' @param Point0 Stima iniziale del centro (vettore lunghezza 2)
#' @param Rad0 Raggio del cerchio (fisso, non stimato)
#' @return Lista con radius, point, mad, ArcCov, conv, rel
#' @export
least_squares_circle_centre <- function(P, Point0, Rad0) {
  par <- c(Point0, Rad0)
  maxiter <- 200L
  iter <- 0L
  conv <- FALSE
  rel  <- TRUE
  
  while (iter < maxiter && !conv && rel) {
    fg   <- func_grad_circle_centre(P, par)
    dist <- fg$dist
    J    <- fg$J
    
    SS0 <- sqrt(sum(dist^2))
    A <- t(J) %*% J
    b <- t(J) %*% dist
    p <- tryCatch(-solve(A, b), error = function(e) NULL)
    if (is.null(p)) {
      rel <- FALSE
      break
    }
    
    par[1:2] <- par[1:2] + as.vector(p)
    
    dist <- func_grad_circle_centre(P, par)$dist
    SS1 <- sqrt(sum(dist^2))
    if (SS1 > SS0) {
      par[1:2] <- par[1:2] - 0.95 * as.vector(p)
      dist <- func_grad_circle_centre(P, par)$dist
      SS1 <- sqrt(sum(dist^2))
    }
    
    if (rcond(A) < 10000 * .Machine$double.eps) {
      rel <- FALSE
    }
    
    if (abs(SS0 - SS1) < 1e-5) {
      conv <- TRUE
    }
    
    iter <- iter + 1L
  }
  
  Point <- par[1:2]
  if (conv && rel) {
    U <- P[, 1] - par[1]
    V <- P[, 2] - par[2]
    ang <- atan2(V, U) + pi
    I <- rep(FALSE, 100)
    ang_idx <- ceiling(ang / (2 * pi) * 100)
    # NOTA FEDELTA' (vedi sopra): nessun clamp per ang_idx==0, fedele
    # all'originale - in R questo semplicemente ignora l'elemento
    # invece di dare errore come farebbe MATLAB.
    valid_idx <- ang_idx >= 1 & ang_idx <= 100
    I[ang_idx[valid_idx]] <- TRUE
    ArcCov <- sum(I) / 100
    d <- sqrt(U^2 + V^2) - Rad0
    mad <- mean(abs(d))
  } else {
    mad <- 0
    ArcCov <- 0
  }
  
  list(
    radius = Rad0,
    point  = Point,
    mad    = mad,
    ArcCov = ArcCov,
    conv   = conv,
    rel    = rel
  )
}
