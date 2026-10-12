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

# VARIABILITA' DEI RISULTATI CON LA BASE ORTONORMALE CASUALE (come nell'originale, che usa rand(3,1) in orthonormal_vectors.m)
# Ripete K volte la corsa di alcuni alberi (seed 0, 8 modelli per albero) con realizzazioni diverse del vettore di partenza pseudo-casuale, che NON tocca il generatore di numeri casuali di R,
# e riassume, per albero, la variazione tra realizzazioni della mediana dei modelli: deviazione standard e coefficiente di variazione (CV) di DBHcyl ufficiale, volume 1.3-4.0 m, 0-4 m, 0-8 m e fusto intero.
# PREVISIONI SCRITTE PRIMA (dopo la prima misura con una realizzazione)
#  P_V1: il CV del volume 1.3-4.0 m (mediana dei modelli di un albero) e' inferiore all'1% nella maggior parte degli alberi (almeno 4 su 6).
#  P_V2: il CV del volume del fusto intero e' maggiore di quello del volume 1.3-4.0 m in tutti gli alberi, con mediana tra gli alberi sopra l'1%.
#  P_V3: il CV del DBHcyl ufficiale e' inferiore allo 0.5% in tutti gli alberi.
# USO (sessione pulita, dopo carica_libreria_pulita(); circa 30 minuti con K = 5):  source("variabilita_base_casuale.R");  v <- variabilita_base_casuale()
variabilita_base_casuale <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", K = 5, alberi = c(4, 16, 35), nuvole = c("alberi_corretti", "alberi_corretti_las11")) {
  stopifnot(exists("raccolta_modelli_quote"), exists("orthonormal_vectors"), exists("cross_product"))
  orig <- get("orthonormal_vectors", envir = globalenv()); on.exit(assign("orthonormal_vectors", orig, envir = globalenv()), add = TRUE)
  colonne <- c("dbh_matlab_cm", "v_13_40", "v_0_4", "v_0_8", "v_stem"); tutti <- list()
  for (k in seq_len(K)) {
    contatore <- 100000L * k
    casuale <- function(U) {
      contatore <<- contatore + 1L; j <- contatore
      V0 <- c(sin(12.9898 * j), sin(78.233 * j + 1), sin(37.719 * j + 2))
      V <- cross_product(V0, U); nV <- sqrt(sum(V^2))
      if (!is.finite(nV) || nV < 1e-10) { V <- cross_product(c(1, 0, 0), U); nV <- sqrt(sum(V^2)) }
      W <- cross_product(V, U); W <- W / sqrt(sum(W^2)); V <- V / nV
      list(V = V, W = W)
    }
    assign("orthonormal_vectors", casuale, envir = globalenv())
    for (cl in nuvole) {
      nuovo <- raccolta_modelli_quote(file.path(base, cl), cartella_R = base, file_out = tempfile(fileext = ".rds"), offsets = 0, alberi = alberi)
      if (is.null(nuovo$modelli) || nrow(nuovo$modelli) == 0) stop("nessun modello prodotto su ", cl, call. = FALSE)
      m <- nuovo$modelli; m$nuvola <- cl; m$realizzazione <- k; tutti[[length(tutti) + 1]] <- m
    }
    cat(sprintf("realizzazione %d di %d fatta\n", k, K))
  }
  d <- do.call(rbind, lapply(tutti, function(m) m[, c("nuvola", "realizzazione", "N", "modello", colonne)]))
  # mediana dei modelli per albero e realizzazione, poi CV tra realizzazioni
  med <- aggregate(d[, colonne], by = list(nuvola = d$nuvola, albero = d$N, realizzazione = d$realizzazione), FUN = median, na.rm = TRUE)
  cv <- aggregate(med[, colonne], by = list(nuvola = med$nuvola, albero = med$albero), FUN = function(x) 100 * sd(x) / abs(mean(x)))
  cat("\nCV (%) tra realizzazioni della mediana dei modelli, per albero:\n"); print(cv, row.names = FALSE, digits = 3)
  cat("\nmediana dei CV tra gli alberi (%):\n"); print(round(sapply(cv[, colonne], median), 3))
  p1 <- mean(cv$v_13_40 < 1) >= 4 / 6; p2 <- all(cv$v_stem > cv$v_13_40) && median(cv$v_stem) > 1; p3 <- all(cv$dbh_matlab_cm < 0.5)
  cat(sprintf("\nP_V1 (CV del volume 1.3-4.0 m sotto l'1%% in almeno 4 alberi su 6): %s\nP_V2 (CV del fusto intero maggiore in tutti gli alberi e mediana sopra l'1%%): %s\nP_V3 (CV del DBHcyl ufficiale sotto lo 0.5%% in tutti): %s\n",
              if (p1) "CONFERMATA" else "NON CONFERMATA", if (p2) "CONFERMATA" else "NON CONFERMATA", if (p3) "CONFERMATA" else "NON CONFERMATA"))
  invisible(list(modelli = d, mediane = med, cv = cv))
}
