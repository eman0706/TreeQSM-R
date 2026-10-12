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

# ============================================================
# PROVA: Gauss-Newton FEDELE al posto di Levenberg-Marquardt in least_squares_cylinder
# Perche': il confronto riga per riga ha mostrato che la least_squares_cylinder di produzione (helpers.R, ripresa identica da
# correzioni_sessione_oggi.R) NON e' quella di TreeQSM 2.4.1. Differenze (originale: least_squares_cylinder.m):
#   1. parametri iniziali 0 (porting: 1e-6);
#   2. aggiornamento SEMPRE accettato, Gauss-Newton puro (porting: smorzamento LM, passo accettato solo se migliora);
#   3. arresto quando |SS0-SS1| < 1e-4 nella STESSA iterazione (porting: 1e-5 contro l'iterazione precedente);
#   4. rcond(-A) < 10000*eps sulla matrice non regolarizzata (porting: A_lm, soglia 1e-12);
#   5. SurfCov calcolata SOLO se rel && conv, altrimenti 0 (porting: calcolata sempre);
#   6. mad = media di |dist| su TUTTI i punti: nell'originale "if nargin >= 6" non e' mai vero (la funzione ha 4 argomenti),
#      quindi il sottoinsieme a peso massimo non viene mai usato (porting: lo usa);
#   7. radius, length, start, axis, mad, SurfCov sono memorizzati in precisione singola (single).
# Il commento di helpers.R dichiarava l'impatto "trascurabile" sulla base di due casi sintetici e dell'albero 47: sull'albero 35 il
# raggio del fusto differisce di 1-1.5 cm e la SurfCov di 0.1-0.2 (originale 0.77-0.90, porting 0.95-0.97).
# PREVISIONI SCRITTE PRIMA DELL'ESECUZIONE (alberi 35, 4, 16, nuvola 15, un seed, 8 modelli ciascuno):
#  P_G1: albero 35: il DBH medio sugli 8 modelli con GN fedele sale di almeno 0.5 cm rispetto alla versione LM (12.0 cm di raggio ->
#        diametro 23.9-24.5 cm) e si avvicina ai 25.3-26.8 cm dell'originale in Octave.
#  P_G2: la differenza mediana per albero tra originale (Octave) e porting scende, sui tre alberi, da 2.1/0.6/3.1 cm (alberi 4/16/35) a
#        meno di 1.0 cm.
#  P_G3: la SurfCov mediana dei cilindri del fusto scende (verso 0.8-0.9).
#  P_G4 (incerta): su tutti i 21 alberi la mediana |errore| resta tra 3.0% e 4.0% (non cambia la conclusione di fondo).
# Se P_G1 e P_G2 non sono confermate, la causa dello scarto NON e' qui e va cercata altrove (parent_cylinder, relative_size,
# surface_coverage_filtering con fit_cerchio_3d, base ortonormale fissa).
# Uso (sessione pulita):
#   source(".../sessione_pulita.R"); carica_libreria_pulita()
#   source(".../prova_gn_fedele.R")
#   r <- prova_gn_fedele(alberi = c(35, 4, 16))                       # alcuni minuti
#   r <- prova_gn_fedele(alberi = NULL)                               # tutti i 21 alberi della nuvola 15, circa 20 minuti
#   (opzionale, per il confronto diretto con l'originale) csv_octave = "F:/.../risultati_octave_alberi_corretti_v2.csv"
# La versione LM viene ripristinata alla fine (on.exit): la prova non modifica nulla su disco.
# ============================================================
a_single <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4), "numeric", size = 4, n = length(x))

least_squares_cylinder_gn <- function(P, cyl0, weight = NULL, Q = NULL) {
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

prova_gn_fedele <- function(cartella_nuvola = "alberi_corretti", base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", alberi = c(35, 4, 16), csv_octave = NULL) {
  stopifnot(exists("raccolta_modelli_quote"), exists("least_squares_cylinder"))
  lm_salvata <- get("least_squares_cylinder", envir = globalenv())
  on.exit(assign("least_squares_cylinder", lm_salvata, envir = globalenv()), add = TRUE)
  assign("least_squares_cylinder", least_squares_cylinder_gn, envir = globalenv())
  cart <- file.path(base, cartella_nuvola); tmp <- tempfile(fileext = ".rds")
  t0 <- Sys.time()
  nu <- raccolta_modelli_quote(cart, cartella_R = base, file_out = tmp, offsets = 0, alberi = alberi)$modelli
  cat(sprintf("\nTempo: %.1f minuti\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rif <- read.csv(file.path(cart, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  lmv <- readRDS(file.path(cart, "modelli_quote_vol.rds"))$modelli; lmv <- lmv[lmv$off == 0, ]
  mm <- merge(nu, lmv, by = c("N", "off", "modello"), suffixes = c("_gn", "_lm"))
  cat("\n=== DBH (cm) per modello: Gauss-Newton fedele (gn) contro Levenberg-Marquardt di produzione (lm) ===\n")
  print(data.frame(N = mm$N, modello = mm$modello, gn = round(mm$dbh_matlab_cm_gn, 2), lm = round(mm$dbh_matlab_cm_lm, 2), diff = round(mm$dbh_matlab_cm_gn - mm$dbh_matlab_cm_lm, 2)), row.names = FALSE)
  med <- function(d, col) median(d[[col]][is.finite(d[[col]]) & d[[col]] >= 5])
  tb <- do.call(rbind, lapply(sort(unique(nu$N)), function(n) { v <- rif$ref_cm[rif$N == n]
    g <- med(nu[nu$N == n, ], "dbh_matlab_cm"); l <- med(lmv[lmv$N == n, ], "dbh_matlab_cm")
    data.frame(N = n, ref_cm = v, DBH_GN = round(g, 1), DBH_LM = round(l, 1), errGN = round(100 * (g - v) / v, 1), errLM = round(100 * (l - v) / v, 1)) }))
  if (!is.null(csv_octave) && file.exists(csv_octave)) {
    o <- read.csv(csv_octave); tb$DBH_Octave <- sapply(tb$N, function(n) round(median(o$dbh_cyl_cm[o$N == n & is.finite(o$dbh_cyl_cm) & o$dbh_cyl_cm >= 5]), 1))
    tb$errOct <- round(100 * (tb$DBH_Octave - tb$ref_cm) / tb$ref_cm, 1)
    tb$GN_meno_Oct <- round(tb$DBH_GN - tb$DBH_Octave, 1); tb$LM_meno_Oct <- round(tb$DBH_LM - tb$DBH_Octave, 1)
  }
  cat("\n=== Per albero (mediana dei modelli, un seed) ===\n"); print(tb, row.names = FALSE)
  cat(sprintf("\nMediana |errore|: GN fedele %.2f%% | LM di produzione %.2f%% | segno mediano GN %+.1f%% LM %+.1f%%\n",
              median(abs(tb$errGN), na.rm = TRUE), median(abs(tb$errLM), na.rm = TRUE), median(tb$errGN, na.rm = TRUE), median(tb$errLM, na.rm = TRUE)))
  if ("GN_meno_Oct" %in% names(tb)) cat(sprintf("Differenza mediana assoluta dall'originale (Octave): GN %.2f cm | LM %.2f cm (P_G2: soglia 1.0 cm)\n", median(abs(tb$GN_meno_Oct), na.rm = TRUE), median(abs(tb$LM_meno_Oct), na.rm = TRUE)))
  invisible(list(modelli = nu, confronto = mm, per_albero = tb))
}
