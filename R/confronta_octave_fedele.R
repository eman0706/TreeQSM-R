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
# Confronto per albero: TreeQSM originale (Octave) contro porting R fedele (modelli_quote_fedele.rds)
# Uso (R, non serve la libreria pulita):
#   source("confronta_octave_fedele.R")
#   c15 <- "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti"          # cartella con riferimento_alberi.csv e rds
#   c11 <- "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti_las11"
#   r15 <- confronta_octave_fedele(c15, c("D:/.../risultati_octave_15_rip2.csv", "D:/.../risultati_octave_15_rip1.csv"), "nuvola 15")
#   r11 <- confronta_octave_fedele(c11, "D:/.../risultati_octave_11_rip1.csv", "nuvola 11")
# Il PRIMO csv e' quello principale (tutte le sue ripetizioni entrano nella mediana Octave per albero);
# i csv successivi servono solo a misurare la ripetibilita' di Octave (stessa ripetizione 1 in due corse).
# Per ogni albero: DBH cilindri = mediana sui modelli (>= 5 cm); R: mediana dei seed (colonna off).
# ============================================================
confronta_octave_fedele <- function(cartella, csv_octave, etichetta = "", file_rds = "modelli_quote_fedele.rds",
                                    nome_rds_colonna = "dbh_matlab_cm") {
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE)
  rif <- rif[rif$fonte_ref == "campo", ]
  R <- readRDS(file.path(cartella, file_rds))$modelli
  ok <- function(x) is.finite(x) & x >= 5
  med <- function(x) { x <- x[ok(x)]; if (length(x)) median(x) else NA_real_ }
  O <- lapply(csv_octave, read.csv, stringsAsFactors = FALSE)
  O1 <- O[[1]]
  alberi <- sort(intersect(intersect(unique(R$N), unique(O1$N)), rif$N))
  # R: mediana sui modelli per (albero, seed), poi sui seed
  r_seed <- aggregate(R[[nome_rds_colonna]], list(N = R$N, off = R$off), med)
  r_alb  <- aggregate(r_seed$x, list(N = r_seed$N), median, na.rm = TRUE)
  o_rip  <- aggregate(O1$dbh_cyl_cm, list(N = O1$N, rip = O1$rip), med)
  o_alb  <- aggregate(o_rip$x, list(N = o_rip$N), median, na.rm = TRUE)
  tb <- data.frame(N = alberi)
  tb$ref <- rif$ref_cm[match(tb$N, rif$N)]
  tb$R   <- round(r_alb$x[match(tb$N, r_alb$N)], 2)
  tb$Oct <- round(o_alb$x[match(tb$N, o_alb$N)], 2)
  tb$R_meno_Oct <- round(tb$R - tb$Oct, 2)
  tb$errR   <- round(100 * (tb$R - tb$ref) / tb$ref, 1)
  tb$errOct <- round(100 * (tb$Oct - tb$ref) / tb$ref, 1)
  cat(sprintf("\n===== %s: porting R fedele (%d seed) contro Octave (%d ripetizione/i) =====\n",
              etichetta, length(unique(R$off)), length(unique(O1$rip))))
  print(tb, row.names = FALSE)
  d <- tb$R_meno_Oct
  cat(sprintf("\nAlberi confrontati: %d\n", sum(is.finite(d))))
  cat(sprintf("Differenza R - Octave: mediana %+.2f cm | mediana |diff| %.2f cm | massima |diff| %.2f cm (albero %s)\n",
              median(d, na.rm = TRUE), median(abs(d), na.rm = TRUE), max(abs(d), na.rm = TRUE), tb$N[which.max(abs(d))]))
  cat(sprintf("Alberi con |diff| < 1 cm: %d / %d | < 2 cm: %d / %d\n", sum(abs(d) < 1, na.rm = TRUE), sum(is.finite(d)),
              sum(abs(d) < 2, na.rm = TRUE), sum(is.finite(d))))
  cat(sprintf("Correlazione R-Octave: Pearson %.3f | Spearman %.3f\n", cor(tb$R, tb$Oct, use = "complete.obs"),
              cor(tb$R, tb$Oct, use = "complete.obs", method = "spearman")))
  cat(sprintf("Mediana |errore| rispetto al campo: R %.2f%% | Octave %.2f%% | segno mediano R %+.1f%% Octave %+.1f%%\n",
              median(abs(tb$errR), na.rm = TRUE), median(abs(tb$errOct), na.rm = TRUE),
              median(tb$errR, na.rm = TRUE), median(tb$errOct, na.rm = TRUE)))
  w <- tryCatch(suppressWarnings(wilcox.test(abs(tb$errR), abs(tb$errOct), paired = TRUE)$p.value), error = function(e) NA)
  cat(sprintf("Wilcoxon appaiato su |errore| R contro Octave: p = %.3f\n", w))
  # Rumore intrinseco: variabilita' tra ripetizioni di Octave e tra seed di R, per albero
  if (length(unique(o_rip$rip)) > 1) {
    w2 <- reshape(o_rip, idvar = "N", timevar = "rip", direction = "wide")
    cc <- grep("^x\\.", names(w2)); dd <- abs(w2[[cc[1]]] - w2[[cc[2]]])
    cat(sprintf("Ripetibilita' Octave (rip %s contro rip %s): mediana |diff| %.2f cm, massima %.2f cm\n",
                sub("x\\.", "", names(w2)[cc[1]]), sub("x\\.", "", names(w2)[cc[2]]), median(dd, na.rm = TRUE), max(dd, na.rm = TRUE)))
  }
  if (length(unique(r_seed$off)) > 1) {
    w3 <- reshape(r_seed, idvar = "N", timevar = "off", direction = "wide")
    cc <- grep("^x\\.", names(w3)); rg <- apply(w3[, cc, drop = FALSE], 1, function(z) diff(range(z, na.rm = TRUE)))
    cat(sprintf("Variabilita' tra seed di R (max-min per albero): mediana %.2f cm, massima %.2f cm\n", median(rg, na.rm = TRUE), max(rg, na.rm = TRUE)))
  }
  if (length(O) > 1) for (k in 2:length(O)) {
    Ok <- O[[k]]; ok_rip <- aggregate(Ok$dbh_cyl_cm, list(N = Ok$N, rip = Ok$rip), med)
    a <- o_rip[o_rip$rip == 1, ]; b <- ok_rip[ok_rip$rip == 1, ]; m <- merge(a, b, by = "N")
    if (nrow(m)) cat(sprintf("Octave rip 1: stessa ripetizione in due corse (file 1 contro file %d): mediana |diff| %.2f cm, massima %.2f cm su %d alberi\n",
                             k, median(abs(m$x.x - m$x.y), na.rm = TRUE), max(abs(m$x.x - m$x.y), na.rm = TRUE), nrow(m)))
  }
  # Livello dei singoli modelli (8 per albero) dove l'ordine e' confrontabile
  invisible(list(per_albero = tb, r_seed = r_seed, o_rip = o_rip))
}
