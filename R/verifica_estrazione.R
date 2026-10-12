# ============================================================
# Verifica della versione R dell'estrazione contro la versione Python (LAS 15, 21 alberi di campo).
# Richiede: dbh_diretto.R, estrazione_alberi.R, i file tronchi_las15_riferimento.csv e posizioni_alberi_las15.csv
# (e riferimento_alberi.csv della v2) nella cartella `cartella`.
# Uso:
#   source(".../dbh_diretto.R"); source(".../estrazione_alberi.R"); source(".../verifica_estrazione.R")
#   res <- verifica_estrazione("F:/.../Data/Area_2_LAS_15.las", "F:/.../R/alberi_corretti")
# Tempo atteso: alcuni minuti.
# ============================================================
verifica_estrazione <- function(percorso_las, cartella) {
  stopifnot(exists("dbh_diretto"), exists("estrai_albero"))
  t0 <- Sys.time()
  las <- carica_las_xyz(percorso_las)
  cat("LAS letto:", nrow(las), "punti\n")
  zrel <- dtm_celle(las$X, las$Y, las$Z)
  tronchi <- rileva_tronchi(las$X, las$Y, zrel)
  ref_t <- read.csv(file.path(cartella, "tronchi_las15_riferimento.csv"))
  cat("Tronchi rilevati in R:", nrow(tronchi), "| in Python:", nrow(ref_t), "\n")
  dmin <- sapply(seq_len(nrow(tronchi)), function(i) min(sqrt((ref_t$x - tronchi$x[i])^2 + (ref_t$y - tronchi$y[i])^2)))
  cat(sprintf("Distanza di ciascun tronco R dal tronco Python piu' vicino: mediana %.3f m, massima %.2f m; entro 0.1 m: %d su %d\n",
              median(dmin), max(dmin), sum(dmin < 0.1), length(dmin)))

  pos <- read.csv(file.path(cartella, "posizioni_alberi_las15.csv"))
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE)
  rif <- rif[rif$fonte_ref == "campo", ]
  pasi15 <- c("1"=1.3,"3"=5.6,"4"=1.1,"5"=3.5,"8"=3.6,"10"=2.9,"14"=2.9,"16"=2.4,"17"=3.3,"20"=4.1,"21"=3.4,"23"=2.8,
              "25"=1.5,"28"=3.0,"30"=0.4,"32"=1.8,"35"=3.0,"39"=3.8,"40"=3.2,"42"=5.1,"45"=0.9)
  atteso_python <- c("1"=22.7,"3"=28.3,"4"=37.7,"5"=35.9,"8"=27.5,"10"=25.6,"14"=33.0,"16"=20.2,"17"=27.8,"20"=24.9,"21"=34.7,
                     "23"=25.2,"25"=31.8,"28"=36.3,"30"=27.5,"32"=48.1,"35"=24.9,"39"=24.6,"40"=30.5,"42"=25.8,"45"=27.1)
  res <- do.call(rbind, lapply(seq_len(nrow(rif)), function(i) {
    n <- rif$N[i]; v <- rif$ref_cm[i]
    p <- pos[pos$N == n, ]
    d <- sqrt((tronchi$x - p$x)^2 + (tronchi$y - p$y)^2); k <- which.min(d)
    if (d[k] > 0.5) return(data.frame(N = n, ref_cm = v, dist_tronco_m = round(d[k], 2), n_punti = NA, n_punti_python = rif$n_punti[i],
                                     dbh_cm = NA, err_segno = NA, originale = pasi15[[as.character(n)]], atteso_python = atteso_python[[as.character(n)]], diff_da_python = NA))
    c0 <- c(tronchi$x[k], tronchi$y[k]); R <- raggio_estrazione(c0, tronchi)
    P <- estrai_albero(las$X, las$Y, zrel, c0, R)
    r <- dbh_diretto(P)
    data.frame(N = n, ref_cm = v, dist_tronco_m = round(d[k], 2), n_punti = nrow(P), n_punti_python = rif$n_punti[i],
               dbh_cm = round(r$diametro_cm, 1), err_segno = round(100 * (r$diametro_cm - v) / v, 1), originale = pasi15[[as.character(n)]],
               atteso_python = atteso_python[[as.character(n)]], diff_da_python = round(r$diametro_cm - atteso_python[[as.character(n)]], 1))
  }))
  print(res, row.names = FALSE)
  cat(sprintf("\nMediana dell'errore assoluto: R %.2f%% | originale MATLAB %.2f%% | media %.2f%% | segno mediano %.1f%% | oltre 15%%: %d | non misurati: %d\n",
              median(abs(res$err_segno), na.rm = TRUE), median(res$originale), mean(abs(res$err_segno), na.rm = TRUE),
              median(res$err_segno, na.rm = TRUE), sum(abs(res$err_segno) > 15, na.rm = TRUE), sum(is.na(res$dbh_cm))))
  cat(sprintf("Differenza massima dal DBH Python: %.2f cm | rapporto mediano punti R / Python: %.3f\n",
              max(abs(res$diff_da_python), na.rm = TRUE), median(res$n_punti / res$n_punti_python, na.rm = TRUE)))
  cat("Tempo:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minuti\n")
  invisible(res)
}
