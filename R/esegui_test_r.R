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
# ESEGUI I TEST R CHE NON RICHIEDONO TREEQSM (circa un minuto)
#   1. DBH diretto a 1.3 m sulle due nuvole, controllato contro la versione Python (P_D1, P_D2)
#   2. DBH diretto alle tre quote di campo (1.30, 2.00, 4.00 m) sulle due nuvole, controllato contro Python (P_R1, P_R2)
# Uso (sessione R nuova va benissimo, non serve il bootstrap del porting):
#   source("F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/esegui_test_r.R"); out <- esegui_test_r()
# Mandare tutto l'output della console.
# ============================================================
esegui_test_r <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R") {
  source(file.path(base, "dbh_diretto.R")); source(file.path(base, "valida_dbh_diretto.R")); source(file.path(base, "dbh_diretto_quote.R"))
  c15 <- file.path(base, "alberi_corretti"); c11 <- file.path(base, "alberi_corretti_las11")
  cat("########## 1. DBH diretto a 1.3 m ##########\n")
  r15 <- valida_dbh_diretto(c15, "15"); r11 <- valida_dbh_diretto(c11, "11")
  cat("\n########## 2. DBH diretto alle tre quote ##########\n")
  q15 <- valida_dbh_diretto_quote(c15, "15"); q11 <- valida_dbh_diretto_quote(c11, "11")
  pd1 <- max(abs(r15$diff_da_python), na.rm = TRUE) < 0.5
  pd2 <- abs(median(abs(r11$err_segno)) - 2.47) <= 0.1 && abs(mean(abs(r11$err_segno)) - 2.23) <= 0.1 && max(abs(r11$diff_da_python), na.rm = TRUE) < 0.5
  pr1 <- q15$diff_max < 0.5 && q11$diff_max < 0.5
  pr2 <- q15$scarto_stat <= 0.1 && q11$scarto_stat <= 0.1
  cat("\n########## ESITI DELLE PREVISIONI ##########\n")
  cat(sprintf("P_D1 (nuvola 15, 1.3 m: differenza massima da Python < 0.5 cm): %s\n", if (pd1) "CONFERMATA" else "NON CONFERMATA"))
  cat(sprintf("P_D2 (nuvola 11, 1.3 m: mediana 2.47%% e media 2.23%% entro 0.1 punti, differenza massima < 0.5 cm): %s\n", if (pd2) "CONFERMATA" else "NON CONFERMATA"))
  cat(sprintf("P_R1 (tre quote: differenza massima da Python < 0.5 cm su entrambe le nuvole): %s\n", if (pr1) "CONFERMATA" else "NON CONFERMATA"))
  cat(sprintf("P_R2 (tre quote: statistiche entro 0.1 punti da Python su entrambe le nuvole): %s\n", if (pr2) "CONFERMATA" else "NON CONFERMATA"))
  invisible(list(r15 = r15, r11 = r11, q15 = q15, q11 = q11))
}
