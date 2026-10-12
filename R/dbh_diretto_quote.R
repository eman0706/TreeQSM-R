# ============================================================
# DBH DIRETTO ALLE TRE QUOTE DI CAMPO (1.30, 2.00, 4.00 m) SULLE DUE NUVOLE, VERSIONE R
# Il file riferimento_alberi.csv contiene i diametri di campo a 1.30 m (ref_cm: media dei due diametri), a 2.00 m (d200) e a 4.00 m (d400).
# I tronchi sono inclinati e le nuvole estratte sono centrate sull'asse a 1.3 m: a 2 e 4 m il fusto e' spostato dall'origine, quindi ogni fetta
# viene ricentrata PRIMA di chiamare dbh_diretto (mediana, poi quattro iterazioni di Kasa entro 0.45 m), come nella versione Python. Le stime fuori da
# 5-80 cm sono considerate non valide e contate. La ricentratura e' stata introdotta dopo una prima esecuzione a 2 e 4 m che dava valori assurdi
# perche' la fetta non era centrata sul fusto.
# Greco et al. (2023, Tabella 9): RMSE e MAPE mediati sulle quote: nuvola 15 0.88 cm e 2.63%; nuvola 11 0.92 cm e 2.73%.
# PREVISIONI SCRITTE PRIMA DELL'ESECUZIONE IN R (i valori Python sono in `atteso` e `att_stat`)
#  P_R1: la versione R differisce da quella Python di al massimo 0.5 cm in ogni stima (a 1.3 m senza ricentratura la differenza massima era 0.10 cm).
#  P_R2: mediana |errore| e MAPE per quota coincidono con quelle Python entro 0.1 punti.
# Uso:  source(".../dbh_diretto.R"); source(".../dbh_diretto_quote.R")
#   r15 <- valida_dbh_diretto_quote(".../alberi_corretti", "15"); r11 <- valida_dbh_diretto_quote(".../alberi_corretti_las11", "11")
# ============================================================
valida_dbh_diretto_quote <- function(cartella, nuvola = c("15", "11")) {
  nuvola <- match.arg(nuvola); stopifnot(exists("dbh_diretto"))
  atteso15 <- list("1.3" = c("1"=22.8, "3"=28.2, "4"=38.0, "5"=35.9, "8"=27.5, "10"=25.6, "14"=33.0, "16"=20.2, "17"=27.8, "20"=24.8, "21"=34.3, "23"=25.2, "25"=31.9, "28"=36.3, "30"=27.5, "32"=49.5, "35"=24.9, "39"=24.6, "40"=30.6, "42"=25.8, "45"=27.0), "2" = c("1"=20.8, "3"=29.0, "4"=35.0, "5"=33.1, "8"=25.4, "10"=25.0, "14"=30.7, "16"=18.5, "17"=25.4, "20"=25.8, "21"=31.9, "23"=24.2, "25"=30.3, "28"=33.9, "30"=32.1, "32"=48.7, "35"=24.1, "39"=24.7, "40"=28.9, "42"=25.5, "45"=26.0), "4" = c("1"=18.5, "3"=26.4, "4"=32.8, "5"=30.2, "8"=24.1, "10"=23.2, "14"=27.9, "16"=17.1, "17"=24.7, "20"=26.5, "21"=29.4, "23"=23.1, "25"=30.3, "28"=30.7, "30"=24.2, "32"=48.7, "35"=22.2, "39"=21.7, "40"=27.7, "42"=25.0, "45"=24.1))
  atteso11 <- list("1.3" = c("1"=22.2, "3"=28.1, "4"=37.4, "5"=35.1, "8"=27.7, "10"=25.5, "14"=32.9, "16"=19.7, "17"=27.0, "20"=25.3, "21"=35.0, "23"=25.9, "25"=31.2, "28"=36.0, "30"=27.5, "32"=49.1, "35"=25.0, "39"=25.2, "40"=31.1, "42"=26.0), "2" = c("1"=20.1, "3"=29.0, "4"=34.4, "5"=32.6, "8"=25.3, "10"=24.8, "14"=30.1, "16"=18.7, "17"=25.0, "20"=25.7, "21"=32.2, "23"=24.4, "25"=30.2, "28"=33.1, "30"=36.5, "32"=47.1, "35"=24.3, "39"=25.8, "40"=29.3, "42"=25.7), "4" = c("1"=17.7, "3"=25.5, "4"=31.4, "5"=30.2, "8"=23.6, "10"=22.7, "14"=27.5, "16"=16.6, "17"=23.7, "20"=23.2, "21"=28.3, "23"=22.2, "25"=29.7, "28"=29.7, "32"=43.2, "35"=21.5, "39"=22.0, "40"=27.8, "42"=23.3))
  stat15 <- list("1.3" = c(1.63, 2.10), "2" = c(1.15, 1.34), "4" = c(2.18, 3.06))
  stat11 <- list("1.3" = c(2.00, 1.95), "2" = c(1.49, 2.09), "4" = c(2.81, 2.92))
  atteso <- if (nuvola == "15") atteso15 else atteso11
  att_stat <- if (nuvola == "15") stat15 else stat11
  greco <- if (nuvola == "15") c(rmse = 0.88, mape = 2.63) else c(rmse = 0.92, mape = 2.73)
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  centra <- function(P, h, half = 0.15) {
    S <- P[abs(P[, 3] - h) < half, , drop = FALSE]; if (nrow(S) < 100) return(NULL)
    c0 <- c(median(S[, 1]), median(S[, 2]))
    for (i in 1:4) {
      m <- sqrt((S[, 1] - c0[1])^2 + (S[, 2] - c0[2])^2) < 0.45; if (sum(m) < 60) return(NULL)
      A <- cbind(2 * S[m, 1], 2 * S[m, 2], 1); b <- S[m, 1]^2 + S[m, 2]^2
      s <- tryCatch(qr.solve(A, b), error = function(e) NULL); if (is.null(s)) return(NULL); c0 <- s[1:2]
    }
    P[, 1] <- P[, 1] - c0[1]; P[, 2] <- P[, 2] - c0[2]; P
  }
  qt <- list("1.3" = list(h = 1.3, col = "ref_cm"), "2" = list(h = 2.0, col = "d200"), "4" = list(h = 4.0, col = "d400"))
  stime <- lapply(rif$N, function(n) {
    P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", n)))); mode(P) <- "numeric"
    sapply(qt, function(q) {
      Q <- centra(P, q$h)
      d <- if (is.null(Q)) NA_real_ else tryCatch(dbh_diretto(Q, h = q$h)$diametro_cm, error = function(e) NA_real_)
      if (is.finite(d) && d >= 5 && d <= 80) d else NA_real_
    })
  })
  stime <- do.call(rbind, stime); rownames(stime) <- rif$N
  cat(sprintf("\n=== nuvola %s ===\n", nuvola)); rm_ <- numeric(0); mp_ <- numeric(0); diff_max <- 0; scarto_stat <- 0
  for (k in names(qt)) {
    ref <- rif[[qt[[k]]$col]]; d <- stime[, k]; ok <- is.finite(ref) & is.finite(d); er <- 100 * (d[ok] - ref[ok]) / ref[ok]
    rmse <- sqrt(mean((d[ok] - ref[ok])^2)); mape <- mean(abs(er)); rm_ <- c(rm_, rmse); mp_ <- c(mp_, mape)
    a <- atteso[[k]]; dm <- max(abs(d[names(a)] - a), na.rm = TRUE); diff_max <- max(diff_max, dm)
    scarto_stat <- max(scarto_stat, abs(median(abs(er)) - att_stat[[k]][1]), abs(mape - att_stat[[k]][2]))
    cat(sprintf("  %.2f m: n = %2d su %d con riferimento | mediana |err| %.2f%% | media |err| %.2f%% | segno mediano %+.1f%% | RMSE %.2f cm | oltre 15%%: %d | stime non valide: %d | diff max da Python %.2f cm\n",
                qt[[k]]$h, sum(ok), sum(is.finite(ref)), median(abs(er)), mape, median(er), rmse, sum(abs(er) > 15), sum(is.finite(ref) & !is.finite(d)), dm))
  }
  cat(sprintf("  media sulle tre quote: RMSE %.2f cm | MAPE %.2f%%   (Greco et al.: %.2f cm | %.2f%%)\n", mean(rm_), mean(mp_), greco["rmse"], greco["mape"]))
  cat(sprintf("  P_R1: differenza massima R - Python = %.2f cm (soglia 0.5) | P_R2: scarto massimo delle statistiche = %.2f punti (soglia 0.1)\n", diff_max, scarto_stat))
  invisible(list(stime = stime, diff_max = diff_max, scarto_stat = scarto_stat, rmse = mean(rm_), mape = mean(mp_)))
}