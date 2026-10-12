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

# least_squares_cylinder_fedele: traduzione di least_squares_cylinder.m (TreeQSM 2.4.1): Gauss-Newton puro, parametri iniziali 0, arresto |SS0-SS1| < 1e-4 nella stessa
# iterazione, rcond(-A) < 10000*eps, SurfCov calcolata solo se rel && conv, mad su tutti i punti, precisione singola dei campi. Sostituisce la versione con
# Levenberg-Marquardt di helpers.R (vedi prova_gn_fedele.R e prova_base_fedele.R: con la base fedele del secondo passaggio la versione LM puo' perdere il fit del piede).
a_single <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4), "numeric", size = 4, n = length(x))

least_squares_cylinder_fedele <- function(P, cyl0, weight = NULL, Q = NULL) {
  res <- 0.03; maxiter <- 50L; iter <- 0L; conv <- FALSE; rel <- TRUE
  NoWeights <- is.null(weight)
  Rot0 <- rotate_to_z_axis(cyl0$axis)$R
  Pt <- (P - matrix(cyl0$start, nrow(P), 3, byrow = TRUE)) %*% t(Rot0)
  par <- c(0, 0, 0, 0, cyl0$radius)
  dist <- rep(0, nrow(P))
  while (iter < maxiter && !conv && rel) {
    fg <- if (NoWeights) func_grad_cylinder(par, Pt) else func_grad_cylinder(par, Pt, weight)
    d0 <- fg$d; J <- fg$J
    SS0 <- sqrt(sum(d0^2))
    A <- t(J) %*% J; b <- t(J) %*% d0
    p <- tryCatch(suppressWarnings(-solve(A, b)), error = function(e) NULL)
    if (is.null(p) || any(!is.finite(p))) { rel <- FALSE; break }
    par <- par + as.vector(p)
    if (rcond(-A) < 10000 * .Machine$double.eps) rel <- FALSE
    dist <- if (NoWeights) func_grad_cylinder(par, Pt)$d else func_grad_cylinder(par, Pt, weight)$d
    SS1 <- sqrt(sum(dist^2))
    if (abs(SS0 - SS1) < 1e-4) conv <- TRUE
    iter <- iter + 1L
  }
  rot2 <- form_rotation_matrices(par[3:4])
  Axis <- as.numeric(t(Rot0) %*% t(rot2$R) %*% c(0, 0, 1))
  Point <- as.numeric(t(Rot0) %*% c(par[1], par[2], 0) + cyl0$start)
  P_use <- if (!is.null(Q) && nrow(Q) > 5) Q else P
  H <- as.numeric(P_use %*% Axis); hmin <- min(H)
  len <- a_single(abs(max(H) - hmin))
  Point <- Point - (sum(Axis * Point) - hmin) * Axis
  mad_val <- mean(abs(dist))   # originale: sempre tutti i punti (nargin >= 6 non e' mai vero)
  cyl <- list(radius = a_single(par[5]), length = len, start = a_single(Point), axis = a_single(Axis),
              mad = a_single(if (is.finite(mad_val)) mad_val else 0), conv = conv, rel = rel, SurfCov = 0)
  if (!any(is.na(Axis)) && !any(is.na(Point)) && rel && conv) {
    nl <- max(3L, ceiling(cyl$length / res)); ns <- min(36L, max(8L, ceiling(2 * pi * cyl$radius / res)))
    sc <- tryCatch(surface_coverage(P_use, Axis, Point, nl, ns, 0.8 * cyl$radius), error = function(e) 0)
    cyl$SurfCov <- a_single(if (is.list(sc)) sc[[1]] else sc)
  }
  cyl
}

