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
# Validazione del DBH diretto sulle nuvole dell'estrazione v2.
# Non serve treeqsm: dura pochi secondi.
# Uso (due passi, nessuna variabile globale da impostare):
#   source("F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/dbh_diretto.R")
#   source("F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/valida_dbh_diretto.R")
#   res <- valida_dbh_diretto("F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti", "15")
#   res11 <- valida_dbh_diretto("F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti_las11", "11")
# Previsioni scritte prima (esecuzione in R): P_D1 nuvola 15: differenza massima dalla versione Python < 0.5 cm; P_D2 nuvola 11: mediana |errore| 2.47% e
# media 2.23% (valori Python) entro 0.1 punti e differenza massima dalla versione Python < 0.5 cm.
# ============================================================
valida_dbh_diretto <- function(cartella, nuvola = c("15", "11")) {
  nuvola <- match.arg(nuvola)
  stopifnot(exists("dbh_diretto"))
  f_rif <- file.path(cartella, "riferimento_alberi.csv")
  cat("Leggo:", f_rif, "| esiste:", file.exists(f_rif), "\n")
  rif <- read.csv(f_rif, stringsAsFactors = FALSE)
  if (!"z_fusto_tracciato" %in% names(rif)) cat("ATTENZIONE: la cartella sembra contenere l'estrazione v1, non la v2\n")
  rif <- rif[rif$fonte_ref == "campo", ]
  cat("Alberi con DBH di campo:", nrow(rif), "\n\n")

  originale15 <- c("1"=1.3,"3"=5.6,"4"=1.1,"5"=3.5,"8"=3.6,"10"=2.9,"14"=2.9,"16"=2.4,"17"=3.3,"20"=4.1,"21"=3.4,"23"=2.8,
              "25"=1.5,"28"=3.0,"30"=0.4,"32"=1.8,"35"=3.0,"39"=3.8,"40"=3.2,"42"=5.1,"45"=0.9)
  # valori ottenuti in Python con lo stesso algoritmo (Hough + LSR entro +-4 cm), per controllare che questa versione R coincida
  atteso15 <- c("1"=22.7,"3"=28.3,"4"=37.7,"5"=35.9,"8"=27.5,"10"=25.6,"14"=33.0,"16"=20.2,"17"=27.8,"20"=24.9,"21"=34.7,
                     "23"=25.2,"25"=31.8,"28"=36.3,"30"=27.5,"32"=48.1,"35"=24.9,"39"=24.6,"40"=30.5,"42"=25.8,"45"=27.1)

  atteso11 <- c("1"=22.3, "3"=27.9, "4"=37.2, "5"=35.1, "8"=27.8, "10"=25.5, "14"=32.9, "16"=19.7, "17"=27.0, "20"=25.2, "21"=34.9, "23"=25.9, "25"=31.2, "28"=36.0, "30"=27.5, "32"=46.9, "35"=25.0, "39"=25.1, "40"=31.1, "42"=26.0)
  # Per la nuvola 11 non ho a disposizione i valori per albero della procedura di Greco et al. nel formato di questo script: la colonna `originale` resta NA.
  originale <- if (nuvola == "15") originale15 else NULL
  atteso <- if (nuvola == "15") atteso15 else atteso11

res <- do.call(rbind, lapply(seq_len(nrow(rif)), function(i) {
    n <- rif$N[i]; v <- rif$ref_cm[i]
    P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", n)))); mode(P) <- "numeric"
    r <- dbh_diretto(P)
    data.frame(N = n, ref_cm = v, dbh_cm = round(r$diametro_cm, 1), kasa_cm = round(r$kasa_cm, 1), hough_cm = round(r$hough_cm, 1),
               err_segno = round(100 * (r$diametro_cm - v) / v, 1),
               originale = if (is.null(originale)) NA_real_ else originale[[as.character(n)]],
               atteso_python = atteso[[as.character(n)]],
               diff_da_python = round(r$diametro_cm - atteso[[as.character(n)]], 1))
  }))
  print(res, row.names = FALSE)
  cat(sprintf("\nMediana dell'errore assoluto: DBH diretto %.2f%%", median(abs(res$err_segno))))
  if (any(is.finite(res$originale))) cat(sprintf(" | originale (procedura di Greco et al.) %.2f%%", median(res$originale)))
  cat(sprintf("\nMedia: %.2f%% | errore con segno mediano %.1f%% | alberi oltre 15%%: %d | NA: %d\n",
              mean(abs(res$err_segno)), median(res$err_segno), sum(abs(res$err_segno) > 15, na.rm = TRUE), sum(is.na(res$dbh_cm))))
  dmax <- max(abs(res$diff_da_python), na.rm = TRUE)
  cat(sprintf("Differenza massima rispetto alla versione Python: %.2f cm (attesa: sotto 0.5 cm)\n", dmax))
  if (nuvola == "11") cat(sprintf("P_D2: scarto della mediana da 2.47%%: %.2f punti, della media da 2.23%%: %.2f punti (soglia 0.1)\n",
                                  abs(median(abs(res$err_segno)) - 2.47), abs(mean(abs(res$err_segno)) - 2.23)))
  invisible(res)
}
