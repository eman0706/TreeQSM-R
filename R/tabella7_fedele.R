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
# TABELLE 6 E 7 CON IL PORTING FEDELE (modelli_quote_fedele.rds, tre seed), alle tre quote di campo (1.30, 2.00, 4.00 m)
# Confronta, sugli stessi alberi (quelli con tutte e tre le stime), il QSM (analogo di DBHcyl alla quota vera dal suolo), il fit diretto
# (valida_dbh_diretto_quote, deterministico, secondi di calcolo) e la procedura di Greco et al. (DiamCyl130/200/400_cm di volumi_autori_las*.csv).
# Per ogni quota e ogni coppia: mediana della differenza tra gli |errori %| (punti; negativo = il primo metodo e' piu' accurato), intervallo di
# fiducia 95% bootstrap sugli alberi (10.000 ricampionamenti, set.seed(1)) e p del Wilcoxon appaiato (esatto = FALSE).
# AGGIUNTE (v2 dello script):
#  - correzione di Holm sui 18 test per nuvola e quota (tabella7 ha la colonna p_holm_18);
#  - ANALISI COMBINATA delle due nuvole: le nuvole sono due ricostruzioni degli STESSI alberi, quindi non sono campioni indipendenti; per ogni albero presente in
#    entrambe si fa la media dell'|errore %| tra le due nuvole e si ripetono i confronti (3 quote x 3 coppie = 9 test, Holm su 9: colonna p_holm_9);
#  - sensibilita' all'esclusione dell'albero 30 (diametro di campo a 2 m maggiore di quello a 1.3 m) alle quote 2 e 4 m, per nuvola.
# Scrive tabella6_fedele.csv, tabella7_fedele.csv, tabella7_combinata_fedele.csv, tabella7_senza_albero30.csv e alberi_errori_fedele.csv nella cartella R.
# Uso (sessione qualsiasi, non serve la libreria TreeQSM):
#   R <- "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R"
#   source(file.path(R, "dbh_diretto.R")); source(file.path(R, "dbh_diretto_quote.R")); source(file.path(R, "tabella7_fedele.R"))
#   t7 <- tabella7_fedele(R)
# ============================================================
tabella7_fedele <- function(R = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", B = 10000, file_rds = "modelli_quote_fedele.rds") {
  stopifnot(exists("valida_dbh_diretto_quote"), exists("dbh_diretto"))
  set.seed(1)
  cfg <- list("15" = list(cart = file.path(R, "alberi_corretti"), au = "volumi_autori_las15.csv"),
              "11" = list(cart = file.path(R, "alberi_corretti_las11"), au = "volumi_autori_las11.csv"))
  hs <- list("1.30" = list(k = "1.3", q = "q130_cyl", ref = "ref_cm", g = "DiamCyl130_cm"),
             "2.00" = list(k = "2", q = "q200_cyl", ref = "d200", g = "DiamCyl200_cm"),
             "4.00" = list(k = "4", q = "q400_cyl", ref = "d400", g = "DiamCyl400_cm"))
  pair <- function(nuvola, quota, nome, a, b) {
    dif <- abs(a) - abs(b)
    bs <- replicate(B, median(dif[sample.int(length(dif), replace = TRUE)]))
    p <- suppressWarnings(wilcox.test(abs(a), abs(b), paired = TRUE, exact = FALSE)$p.value)
    data.frame(nuvola = nuvola, quota = quota, n = length(dif), coppia = nome, diff_med = median(dif), ic_inf = unname(quantile(bs, .025)), ic_sup = unname(quantile(bs, .975)), p = p, stringsAsFactors = FALSE)
  }
  tre_coppie <- function(nuvola, quota, eq, ed, eg) rbind(pair(nuvola, quota, "QSM-Diretto", eq, ed), pair(nuvola, quota, "QSM-Greco", eq, eg), pair(nuvola, quota, "Diretto-Greco", ed, eg))
  r6 <- list(); r7 <- list(); alb <- list()
  for (nv in names(cfg)) {
    cart <- cfg[[nv]]$cart
    rif <- read.csv(file.path(cart, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
    fed <- readRDS(file.path(cart, file_rds))$modelli
    au <- read.csv(file.path(cart, cfg[[nv]]$au), stringsAsFactors = FALSE)
    stopifnot(all(c("DiamCyl130_cm", "DiamCyl200_cm", "DiamCyl400_cm") %in% names(au)))
    stime_dir <- valida_dbh_diretto_quote(cart, nv)$stime     # stampa le statistiche del diretto: devono coincidere con la Tabella 6 della v6
    pa <- function(col, n) { d <- fed[fed$N == n & is.finite(fed[[col]]) & fed[[col]] >= 5, ]; if (nrow(d) == 0) return(NA_real_); median(tapply(d[[col]], d$off, median)) }
    mape_h <- list(Diretto = numeric(0), QSM = numeric(0), Greco = numeric(0)); rmse_h <- mape_h
    cat(sprintf("\n=================== NUVOLA %s ===================\n", nv))
    for (hn in names(hs)) {
      h <- hs[[hn]]
      ref <- rif[[h$ref]]
      q <- sapply(rif$N, function(n) pa(h$q, n))
      d <- stime_dir[as.character(rif$N), h$k]
      g <- au[[h$g]][match(rif$N, au$TreeID)]
      ok <- is.finite(ref) & is.finite(q) & is.finite(d) & is.finite(g)
      er <- list(Diretto = 100 * (d[ok] - ref[ok]) / ref[ok], QSM = 100 * (q[ok] - ref[ok]) / ref[ok], Greco = 100 * (g[ok] - ref[ok]) / ref[ok])
      val <- list(Diretto = d[ok], QSM = q[ok], Greco = g[ok])
      for (m in names(er)) {
        rm <- sqrt(mean((val[[m]] - ref[ok])^2)); mp <- mean(abs(er[[m]]))
        mape_h[[m]] <- c(mape_h[[m]], mp); rmse_h[[m]] <- c(rmse_h[[m]], rm)
        r6[[length(r6) + 1]] <- data.frame(nuvola = nv, metodo = m, quota = hn, n = sum(ok), med_abs = median(abs(er[[m]])), media_abs = mp, segno_med = median(er[[m]]), rmse_cm = rm, stringsAsFactors = FALSE)
      }
      r7[[length(r7) + 1]] <- tre_coppie(nv, hn, er$QSM, er$Diretto, er$Greco)
      alb[[length(alb) + 1]] <- data.frame(nuvola = nv, quota = hn, N = rif$N[ok], e_QSM = er$QSM, e_Diretto = er$Diretto, e_Greco = er$Greco, stringsAsFactors = FALSE)
    }
    for (m in names(mape_h)) r6[[length(r6) + 1]] <- data.frame(nuvola = nv, metodo = m, quota = "media quote", n = NA, med_abs = NA, media_abs = mean(mape_h[[m]]), segno_med = NA, rmse_cm = mean(rmse_h[[m]]), stringsAsFactors = FALSE)
  }
  t6 <- do.call(rbind, r6); t7 <- do.call(rbind, r7); AL <- do.call(rbind, alb)
  t7$p_holm_18 <- p.adjust(t7$p, "holm")
  t7$testo <- sprintf("%+.2f [%+.2f, %+.2f]; p = %.3f (Holm %.3f)", t7$diff_med, t7$ic_inf, t7$ic_sup, t7$p, t7$p_holm_18)
  cat("\n\n############ TABELLA 6 (per metodo, stessi alberi) ############\n"); print(transform(t6, med_abs = round(med_abs, 2), media_abs = round(media_abs, 2), segno_med = round(segno_med, 1), rmse_cm = round(rmse_cm, 2)), row.names = FALSE)
  cat("\n############ TABELLA 7 (differenza mediana tra |errori| in punti, IC95 bootstrap, p Wilcoxon e p corretto con Holm su 18 test) ############\n")
  print(t7[, c("nuvola", "quota", "n", "coppia", "testo")], row.names = FALSE)
  
  # ---------- analisi combinata delle due nuvole (media per albero tra le nuvole) ----------
  comb <- list()
  for (hn in names(hs)) {
    a <- AL[AL$quota == hn & AL$nuvola == "15", ]; b <- AL[AL$quota == hn & AL$nuvola == "11", ]
    com <- intersect(a$N, b$N); a <- a[match(com, a$N), ]; b <- b[match(com, b$N), ]
    mm <- function(x, y) (abs(x) + abs(y)) / 2      # |errore| medio tra le due nuvole (positivo)
    eq <- mm(a$e_QSM, b$e_QSM); ed <- mm(a$e_Diretto, b$e_Diretto); eg <- mm(a$e_Greco, b$e_Greco)
    comb[[length(comb) + 1]] <- tre_coppie("combinata", hn, eq, ed, eg)
  }
  tc <- do.call(rbind, comb); tc$p_holm_9 <- p.adjust(tc$p, "holm")
  tc$testo <- sprintf("%+.2f [%+.2f, %+.2f]; p = %.3f (Holm %.3f)", tc$diff_med, tc$ic_inf, tc$ic_sup, tc$p, tc$p_holm_9)
  cat("\n############ ANALISI COMBINATA (media per albero dell'|errore| tra le due nuvole; Holm su 9 test) ############\n")
  print(tc[, c("quota", "n", "coppia", "testo")], row.names = FALSE)
  
  # ---------- sensibilita' all'albero 30 (quote 2 e 4 m) ----------
  ss <- list()
  for (nv in names(cfg)) for (hn in c("2.00", "4.00")) {
    x <- AL[AL$nuvola == nv & AL$quota == hn & AL$N != 30, ]
    if (nrow(x) >= 6) ss[[length(ss) + 1]] <- tre_coppie(nv, hn, x$e_QSM, x$e_Diretto, x$e_Greco)
  }
  ts <- do.call(rbind, ss); ts$testo <- sprintf("%+.2f [%+.2f, %+.2f]; p = %.3f", ts$diff_med, ts$ic_inf, ts$ic_sup, ts$p)
  cat("\n############ SENZA L'ALBERO 30 (quote 2 e 4 m; non corretto per test multipli) ############\n")
  print(ts[, c("nuvola", "quota", "n", "coppia", "testo")], row.names = FALSE)
  
  write.csv(t6, file.path(R, "tabella6_fedele.csv"), row.names = FALSE); write.csv(t7, file.path(R, "tabella7_fedele.csv"), row.names = FALSE)
  write.csv(tc, file.path(R, "tabella7_combinata_fedele.csv"), row.names = FALSE); write.csv(ts, file.path(R, "tabella7_senza_albero30.csv"), row.names = FALSE)
  write.csv(AL, file.path(R, "alberi_errori_fedele.csv"), row.names = FALSE)
  cat(sprintf("\nScritti in %s: tabella6_fedele.csv, tabella7_fedele.csv, tabella7_combinata_fedele.csv, tabella7_senza_albero30.csv, alberi_errori_fedele.csv\n", R))
  invisible(list(tabella6 = t6, tabella7 = t7, combinata = tc, senza30 = ts, alberi = AL))
}