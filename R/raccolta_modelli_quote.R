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
# DIAMETRO DEL QSM ALLE TRE QUOTE DI CAMPO (1.30, 2.00, 4.00 m) SULLE DUE NUVOLE
# Riesegue treeqsm (stessa configurazione e stessi seed della validazione finale) e legge, per ogni modello, il diametro del fusto
# alle quote vere dal suolo H = 1.30, 2.00, 4.00 m con l'analogo di DBHcyl: cilindro del fusto che attraversa la quota H, poi nuovo fit ai
# minimi quadrati dei punti del fusto con |z - H| < 0.2 m (stesse regole di accettazione di dbh_cylinder: raggio entro 0.8-1.2 volte quello del
# cilindro, asse allineato > 0.9, fit convergente e affidabile; altrimenti diametro del cilindro). Il file riferimento_alberi.csv contiene i
# diametri di campo: ref_cm (1.30 m, media di due diametri), d200 (2.00 m), d400 (4.00 m).
# ATTENZIONE: e' l'ANALOGO di dbh_cylinder alle quote vere dal suolo, non la funzione di MATLAB: dbh_cylinder legge a 1.3 m di lunghezza di
# fusto dall'inizio della catena (circa 0.10 m sopra il suolo), per questo il suo risultato ufficiale (colonna dbh_matlab_cm, rifatta qui come
# controllo) differisce di circa -1.1% ogni 0.1 m dal diametro a 1.30 m vero dal suolo. Il risultato principale del lavoro resta dbh_matlab_cm.
# PREVISIONI SCRITTE PRIMA DELL'ESECUZIONE
#  P_Q1: riproducibilita': la colonna dbh_matlab_cm di questa corsa coincide con quella dei file modelli_v5_riparazione.rds (nuvola 15) e
#        modelli_las11_v5_riparazione.rds (nuvola 11) entro 0.05 cm in almeno il 95% dei modelli (stessi seed, stesso codice).
#  P_Q2: mediana |errore| del diametro del QSM alle quote vere <= 5% a 2.00 m e <= 6% a 4.00 m su entrambe le nuvole.
#  P_Q3: media sulle tre quote di MAPE e di RMSE del QSM peggiore (piu' alta) di quella di Greco et al. (nuvola 15: 0.88 cm e 2.63%; nuvola 11:
#        0.92 cm e 2.73%) su entrambe le nuvole, coerentemente con il risultato a 1.3 m.
#  P_Q4: la MAPE media del diretto (2.17% nuvola 15, 2.32% nuvola 11) resta piu' bassa di quella del QSM.
# VOLUME DEL FUSTO (stessa corsa, nessun tempo in piu'): per ogni modello vengono registrati anche il volume del fusto (catena di cilindri con
# branch == 1) tra 1.3 e 4.0 m, da 0 a 4 m, da 0 a 8 m e totale (litri; i cilindri a cavallo di una quota sono ripartiti in proporzione all'altezza), e i volumi di
# treedata (totale, tronco, rami; unita' come restituite dal porting). Riferimenti: (i) VALIDAZIONE del tratto 1.3-4.0 m contro il volume ricavato dai diametri
# di campo a 1.30, 2.00 e 4.00 m con la formula del tronco di cono (alberi con d400: 14 nuvola 15, 13 nuvola 11); (ii) CONFRONTO TRA PROCEDURE con i volumi della
# procedura di Greco et al. (Volume4m, Volume8m, VolumeStem; file volumi_autori_las15.csv e volumi_autori_las11.csv da copiare nelle cartelle delle nuvole): nessuno dei
# due e' una verita' di campo. NON esistono volumi o altezze misurati in campo: il volume dei rami e il volume totale non sono validati.
# PREVISIONI SUI VOLUMI, scritte prima (08/10/2026; la P_V2 e' stata verificata in anteprima in Python con i valori del DBH diretto: mediana |err| 3.0% nuvola 15, 2.0% nuvola 11)
#  P_V1: volume QSM del tratto 1.3-4.0 m: mediana |errore| <= 10% e segno mediano negativo su entrambe le nuvole (diametri sottostimati di circa il 3.5%).
#  P_V2: volume del tratto 1.3-4.0 m dal DBH diretto: mediana |errore| <= 5% su entrambe le nuvole.
#  P_V3: volume QSM da 0 a 4 m entro +-15% di Volume4m della procedura di Greco et al. in almeno il 75% degli alberi, su entrambe le nuvole.
#  P_V4: la differenza relativa mediana tra QSM e procedura di Greco et al. e' maggiore per il volume da 0 a 8 m che per quello da 0 a 4 m (nuvole piu' rade in alto).
# Tempo: circa 50 minuti (nuvola 15) e 80 minuti (nuvola 11) con offsets = 0:2; controllo veloce: alberi = c(4, 16, 35), offsets = 0 (alcuni minuti).
# Prerequisiti: sessione con la configurazione finale (surface_coverage corretta, define_main_branches fedele attivo, correzioni caricate):
#   source(".../prova_tree_sets_fedele.R"); assign("define_main_branches", define_main_branches_fedele, envir = globalenv())
# Uso:
#   source(".../raccolta_modelli_quote.R")
#   r15 <- raccolta_modelli_quote(c15); r11 <- raccolta_modelli_quote(c11)
#   a15 <- analizza_quote(r15, c15, "F:/.../alberi_corretti/modelli_v5_riparazione.rds", "15")
#   a11 <- analizza_quote(r11, c11, "F:/.../alberi_corretti_las11/modelli_las11_v5_riparazione.rds", "11")
#   # volumi: copiare prima volumi_autori_las15.csv in c15 e volumi_autori_las11.csv in c11
#   v15 <- analizza_volumi(r15, c15, file.path(c15, "volumi_autori_las15.csv"), "15")
#   v11 <- analizza_volumi(r11, c11, file.path(c11, "volumi_autori_las11.csv"), "11")
# ============================================================
dbh_cyl_quota <- function(m, H, semi = 0.2) {
  na <- c(qsm = NA_real_, cyl = NA_real_)
  cy <- m$cylinder; st <- which(cy$branch == 1); if (length(st) == 0) return(na)
  z0 <- cy$start[st, 3]; z1 <- z0 + cy$axis[st, 3] * cy$length[st]
  j <- which(pmin(z0, z1) <= H & pmax(z0, z1) >= H); if (length(j) == 0) return(na)
  i <- st[j[1]]; dq <- 2 * cy$radius[i]; dc <- dq; tp <- m$trunk_points
  if (!is.null(tp) && nrow(tp) > 0) {
    I <- tp[, 3] > H - semi & tp[, 3] < H + semi
    if (sum(I) > 100) {
      cyl0 <- select_cylinders(cy, i)
      cyl <- tryCatch(least_squares_cylinder(tp[I, , drop = FALSE], cyl0), error = function(e) NULL)
      if (!is.null(cyl)) {
        ok <- (2 * cyl$radius > 0.8 * dq) && (2 * cyl$radius < 1.2 * dq); al <- abs(sum(cy$axis[i, ] * cyl$axis))
        if (isTRUE(ok) && is.finite(al) && al > 0.9 && isTRUE(cyl$conv) && isTRUE(cyl$rel)) dc <- 2 * cyl$radius
      }
    }
  }
  c(qsm = 100 * dq, cyl = 100 * dc)
}
vol_fusto <- function(m, z1, z2) {
  cy <- m$cylinder; st <- which(cy$branch == 1); if (length(st) == 0) return(NA_real_)
  za <- cy$start[st, 3]; zb <- za + cy$axis[st, 3] * cy$length[st]
  lo <- pmin(za, zb); hi <- pmax(za, zb); v <- pi * cy$radius[st]^2 * cy$length[st]
  span <- pmax(hi - lo, 1e-9); ov <- pmax(0, pmin(hi, z2) - pmax(lo, z1))
  1000 * sum(v * pmin(1, ov / span))
}
raccolta_modelli_quote <- function(cartella, cartella_R = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", file_out = NULL, offsets = 0:2, alberi = NULL) {
  stopifnot(exists("treeqsm"), exists("inputs_prod"), exists("dbh_cylinder"), exists("select_cylinders"), exists("least_squares_cylinder"))
  if (!exists("define_main_branches_fedele") || !identical(get("define_main_branches", envir = globalenv()), define_main_branches_fedele))
    stop("define_main_branches fedele non attiva: source('.../prova_tree_sets_fedele.R') e assign('define_main_branches', define_main_branches_fedele, envir = globalenv())")
  if (!("full" %in% names(formals(surface_coverage)))) stop("surface_coverage non e' la versione corretta")
  if (is.null(file_out)) file_out <- file.path(cartella, "modelli_quote.rds")
  rs_salvata <- get("relative_size", envir = globalenv()); on.exit(assign("relative_size", rs_salvata, envir = globalenv()), add = TRUE)
  source(file.path(cartella_R, "relative_size.R"))
  if (!any(grepl("ChildSegment", deparse(get("relative_size", envir = globalenv()))))) stop("relative_size.R non ha sostituito la versione attiva")
  treeqsm_orig <- get("treeqsm", envir = globalenv()); on.exit(assign("treeqsm", treeqsm_orig, envir = globalenv()), add = TRUE)
  txt <- readLines(file.path(cartella_R, "treeqsm.R"), warn = FALSE); i_td <- grep("treedata = tdata\\$treedata,", txt)
  if (length(i_td) != 1) stop("treeqsm.R: punto di inserimento non trovato (righe: ", length(i_td), ")")
  txt[i_td] <- paste0(txt[i_td], "\n            trunk_points = trunk_points,"); tf <- tempfile(fileext = ".R"); writeLines(txt, tf); source(tf)
  if (!any(grepl("trunk_points = trunk_points", deparse(get("treeqsm", envir = globalenv()))))) stop("patch di treeqsm non applicata")
  assign(".off", 0L, envir = globalenv())
  assign("set.seed", function(seed, ...) base::set.seed(as.integer((as.numeric(seed) + get(".off", envir = globalenv()) * 7919) %% 2147483647), ...), envir = globalenv())
  on.exit(if (exists("set.seed", envir = globalenv(), inherits = FALSE)) rm("set.seed", envir = globalenv()), add = TRUE)
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  if (!is.null(alberi)) rif <- rif[rif$N %in% alberi, ]
  ris <- if (file.exists(file_out)) readRDS(file_out) else list(modelli = NULL)
  fatti <- if (!is.null(ris$modelli)) unique(ris$modelli$N) else integer(0)
  gv <- function(x) if (is.null(x) || length(x) == 0 || !is.finite(x[1])) NA_real_ else as.numeric(x[1])
  dbh_fed <- function(m) {
    if (is.null(m$trunk_points)) return(NA_real_); cy <- m$cylinder
    r <- tryCatch(dbh_cylinder(cy$branch == 1, cy, seq_along(cy$radius), m$trunk_points), error = function(e) NULL)
    if (is.null(r)) NA_real_ else 100 * r$DBHcyl
  }
  for (i in seq_len(nrow(rif))) {
    n <- rif$N[i]; if (n %in% fatti) next
    P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", n)))); mode(P) <- "numeric"
    Pp <- P[P[, 3] >= 0.15, , drop = FALSE]; if (nrow(Pp) < 500) Pp <- P
    for (o in offsets) {
      assign(".off", as.integer(o), envir = globalenv()); q <- NULL
      errmsg <- NULL; invisible(capture.output(q <- tryCatch(treeqsm(Pp, inputs_prod), error = function(e) { errmsg <<- conditionMessage(e); NULL })))
      if (is.null(q) || length(q) == 0) cat(sprintf("  ATTENZIONE: nessun modello per l'albero %d (seed %d): %s\n", n, o, if (is.null(errmsg)) "treeqsm non ha restituito modelli" else errmsg))
      if (!is.null(q) && length(q) > 0) {
        mod <- do.call(rbind, lapply(seq_along(q), function(k) {
          m <- q[[k]]; a <- dbh_cyl_quota(m, 1.3); b <- dbh_cyl_quota(m, 2.0); d <- dbh_cyl_quota(m, 4.0)
          data.frame(N = n, off = o, modello = k, dbh_matlab_cm = dbh_fed(m), q130_qsm = a[["qsm"]], q130_cyl = a[["cyl"]], q200_qsm = b[["qsm"]], q200_cyl = b[["cyl"]],
                     q400_qsm = d[["qsm"]], q400_cyl = d[["cyl"]],
                     v_13_40 = vol_fusto(m, 1.3, 4.0), v_0_4 = vol_fusto(m, -1, 4.0), v_0_8 = vol_fusto(m, -1, 8.0), v_stem = vol_fusto(m, -1, 1e6),
                     vol_tot_td = gv(m$treedata$TotalVolume), vol_trunk_td = gv(m$treedata$TrunkVolume), vol_branch_td = gv(m$treedata$BranchVolume),
                     n_cil = length(m$cylinder$radius), n_cil_tronco = sum(m$cylinder$branch == 1), h_qsm = max(pmax(m$cylinder$start[, 3], m$cylinder$start[, 3] + m$cylinder$axis[, 3] * m$cylinder$length)), z_ini_tronco = min(m$cylinder$start[m$cylinder$branch == 1, 3]), z_max_nuvola = max(Pp[, 3]), stringsAsFactors = FALSE)
        }))
        ris$modelli <- rbind(ris$modelli, mod)
      }
    }
    saveRDS(ris, file_out); cat(sprintf("[%s] albero %d fatto (%d modelli)\n", format(Sys.time(), "%H:%M:%S"), n, if (is.null(ris$modelli)) 0L else sum(ris$modelli$N == n)))
  }
  invisible(ris)
}
analizza_quote <- function(ris, cartella, rds_precedente = NULL, nuvola = c("15", "11")) {
  nuvola <- match.arg(nuvola)
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  greco <- if (nuvola == "15") c(rmse = 0.88, mape = 2.63) else c(rmse = 0.92, mape = 2.73)
  diretto_mape <- if (nuvola == "15") 2.17 else 2.32
  M <- ris$modelli; per_albero <- function(col, n) {
    d <- M[M$N == n & is.finite(M[[col]]) & M[[col]] >= 5, ]; if (nrow(d) == 0) return(NA_real_)
    median(tapply(d[[col]], d$off, median))
  }
  qt <- list("1.3" = c("q130_cyl", "ref_cm"), "2" = c("q200_cyl", "d200"), "4" = c("q400_cyl", "d400"))
  cat(sprintf("\n=== nuvola %s: diametro del QSM (analogo di DBHcyl) alle quote vere dal suolo ===\n", nuvola)); rm_ <- numeric(0); mp_ <- numeric(0); st <- list()
  for (k in names(qt)) {
    d <- sapply(rif$N, function(n) per_albero(qt[[k]][1], n)); ref <- rif[[qt[[k]][2]]]; ok <- is.finite(d) & is.finite(ref); er <- 100 * (d[ok] - ref[ok]) / ref[ok]
    rmse <- sqrt(mean((d[ok] - ref[ok])^2)); mape <- mean(abs(er)); rm_ <- c(rm_, rmse); mp_ <- c(mp_, mape); st[[k]] <- median(abs(er))
    cat(sprintf("  %.2f m: n = %2d su %d con riferimento | mediana |err| %.2f%% | media |err| %.2f%% | segno mediano %+.1f%% | RMSE %.2f cm | oltre 15%%: %d\n",
                as.numeric(k), sum(ok), sum(is.finite(ref)), median(abs(er)), mape, median(er), rmse, sum(abs(er) > 15)))
  }
  cat(sprintf("  media sulle tre quote: RMSE %.2f cm | MAPE %.2f%%   (Greco et al.: %.2f cm | %.2f%%; diretto: MAPE %.2f%%)\n", mean(rm_), mean(mp_), greco["rmse"], greco["mape"], diretto_mape))
  rip <- NA
  if (!is.null(rds_precedente) && file.exists(rds_precedente)) {
    Pm <- readRDS(rds_precedente)$modelli; mm <- merge(M[, c("N", "off", "modello", "dbh_matlab_cm")], Pm[, c("N", "off", "modello", "dbh_matlab_cm")], by = c("N", "off", "modello"), suffixes = c("_nuovo", "_prima"))
    dd <- abs(mm$dbh_matlab_cm_nuovo - mm$dbh_matlab_cm_prima); rip <- mean(dd <= 0.05, na.rm = TRUE)
    cat(sprintf("  P_Q1 riproducibilita': %d modelli confrontati, %.1f%% entro 0.05 cm (soglia 95%%); differenza massima %.2f cm\n", nrow(mm), 100 * rip, max(dd, na.rm = TRUE)))
  }
  cat(sprintf("  P_Q2: mediana |err| a 2.00 m %.2f%% (soglia 5) e a 4.00 m %.2f%% (soglia 6) | P_Q3: MAPE %.2f%% e RMSE %.2f cm contro %.2f%% e %.2f cm di Greco | P_Q4: MAPE QSM %.2f%% contro diretto %.2f%%\n",
              st[["2"]], st[["4"]], mean(mp_), mean(rm_), greco["mape"], greco["rmse"], mean(mp_), diretto_mape))
  invisible(list(rmse = mean(rm_), mape = mean(mp_), riproducibilita = rip))
}

analizza_volumi <- function(ris, cartella, csv_autori, nuvola = c("15", "11")) {
  nuvola <- match.arg(nuvola); M <- ris$modelli
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  au <- read.csv(csv_autori, stringsAsFactors = FALSE)
  per_albero <- function(col, n) { d <- M[M$N == n & is.finite(M[[col]]) & M[[col]] > 0, ]; if (nrow(d) == 0) return(NA_real_); median(tapply(d[[col]], d$off, median)) }
  tc <- function(d1, d2, h) { d1 <- d1 / 100; d2 <- d2 / 100; pi * h / 12 * (d1^2 + d1 * d2 + d2^2) * 1000 }      # litri
  vf <- tc(rif$ref_cm, rif$d200, 0.7) + tc(rif$d200, rif$d400, 2.0)                                              # tratto 1.3-4.0 m dai diametri di campo
  vq <- sapply(rif$N, function(n) per_albero("v_13_40", n)); ok <- is.finite(vf) & is.finite(vq); er <- 100 * (vq[ok] - vf[ok]) / vf[ok]
  cat(sprintf("\n=== nuvola %s: volume del fusto ===\n", nuvola))
  cat(sprintf("  Tratto 1.3-4.0 m, QSM contro il campo (tronco di cono dai diametri di campo): n = %d | mediana |err| %.1f%% | media |err| %.1f%% | segno mediano %+.1f%% | massimo |err| %.1f%%\n",
              sum(ok), median(abs(er)), mean(abs(er)), median(er), max(abs(er))))
  rr <- list(); for (k in list(c("Volume4m_L", "v_0_4", "0-4 m"), c("Volume8m_L", "v_0_8", "0-8 m"), c("VolumeStem_L", "v_stem", "fusto intero"))) {
    a <- au[[k[1]]][match(rif$N, au$TreeID)]; q <- sapply(rif$N, function(n) per_albero(k[2], n)); g <- is.finite(a) & is.finite(q) & a > 0; rel <- 100 * (q[g] - a[g]) / a[g]
    rr[[k[3]]] <- c(n = sum(g), med_abs = median(abs(rel)), segno = median(rel), entro15 = mean(abs(rel) <= 15), rapporto = median(q[g] / a[g]))
    cat(sprintf("  Fusto %-12s QSM contro procedura di Greco et al. (non e' un riferimento di campo): n = %2d | rapporto mediano %.2f | differenza relativa mediana |.| %.1f%% | con segno %+.1f%% | entro +-15%%: %.0f%%\n",
                k[3], sum(g), median(q[g] / a[g]), median(abs(rel)), median(rel), 100 * mean(abs(rel) <= 15)))
  }
  cat(sprintf("  Volumi di treedata (non validati; unita' del porting): mediana totale %.1f | tronco %.1f | rami %.1f\n",
              median(sapply(rif$N, function(n) per_albero("vol_tot_td", n)), na.rm = TRUE), median(sapply(rif$N, function(n) per_albero("vol_trunk_td", n)), na.rm = TRUE),
              median(sapply(rif$N, function(n) per_albero("vol_branch_td", n)), na.rm = TRUE)))
  cat(sprintf("  P_V1: mediana |err| %.1f%% (soglia 10) e segno mediano %+.1f%% (atteso negativo) | P_V3: entro +-15%% sul volume 0-4 m: %.0f%% (soglia 75%%) | P_V4: differenza mediana 0-8 m %.1f%% contro 0-4 m %.1f%%\n",
              median(abs(er)), median(er), 100 * rr[["0-4 m"]][["entro15"]], rr[["0-8 m"]][["med_abs"]], rr[["0-4 m"]][["med_abs"]]))
  invisible(list(err_sezione = er, confronto = rr))
}
tabella_quote <- function(ris, cartella) {
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  M <- ris$modelli
  pa <- function(col, n) { d <- M[M$N == n & is.finite(M[[col]]) & M[[col]] >= 5, ]; if (nrow(d) == 0) return(NA_real_); median(tapply(d[[col]], d$off, median)) }
  data.frame(N = rif$N, ref130 = rif$ref_cm, qsm130 = round(sapply(rif$N, function(n) pa("q130_cyl", n)), 2), ref200 = rif$d200, qsm200 = round(sapply(rif$N, function(n) pa("q200_cyl", n)), 2),
             ref400 = rif$d400, qsm400 = round(sapply(rif$N, function(n) pa("q400_cyl", n)), 2))
}
tabella_volumi <- function(ris, cartella, csv_autori) {
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  au <- read.csv(csv_autori, stringsAsFactors = FALSE); M <- ris$modelli
  pa <- function(col, n) { d <- M[M$N == n & is.finite(M[[col]]) & M[[col]] > 0, ]; if (nrow(d) == 0) return(NA_real_); median(tapply(d[[col]], d$off, median)) }
  tc <- function(d1, d2, h) { d1 <- d1 / 100; d2 <- d2 / 100; pi * h / 12 * (d1^2 + d1 * d2 + d2^2) * 1000 }
  vf <- tc(rif$ref_cm, rif$d200, 0.7) + tc(rif$d200, rif$d400, 2.0)
  v1340 <- sapply(rif$N, function(n) pa("v_13_40", n)); v04 <- sapply(rif$N, function(n) pa("v_0_4", n)); a4 <- au$Volume4m_L[match(rif$N, au$TreeID)]
  data.frame(N = rif$N, vol_campo_1.3_4.0_L = round(vf, 1), vol_qsm_1.3_4.0_L = round(v1340, 1), err_pct = round(100 * (v1340 - vf) / vf, 1),
             vol_qsm_0_4_L = round(v04, 1), autori_Volume4m_L = round(a4, 1), rapporto_qsm_autori = round(v04 / a4, 2))
}