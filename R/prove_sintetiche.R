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
# PROVE SINTETICHE CON IL CODICE CORRETTO (relative_size fedele + surface_coverage corretta)
#
# I generatori usati in passato (gen_tronco, gen2) esistevano solo nello snapshot: questo script li riscrive dalle definizioni
# del manoscritto (cilindri verticali alti 8 m, densita' 2000 punti/m^2 di superficie, rumore radiale, copertura angolare
# parziale, punti interni, guscio esterno fino a 0.15 m dalla superficie tra 0 e 3 m, tronco vicino di raggio 0.15 m a 0.9 m,
# superficie duplicata e spostata). Le frazioni sono rispetto al numero di punti di superficie. Come per gli alberi veri, la
# nuvola parte da z = 0.15 m e le stime sono: A = dbh_robusto (lettura personalizzata), F = DBHcyl fedele di MATLAB (mediana
# dei modelli), E = diametro della catena a 1.3 m senza rifit, Diretto = dbh_diretto sulla nuvola.
#
# CONDIZIONI: 8 ideali (raggio 0.12 e 0.25 m x rumore 0.005 e 0.010 m x copertura 1.0 e 0.6) e 7 scenari su un fusto di 24 cm
# (rumore 0.005, copertura 1): interni 0.2 e 0.5; guscio 0.5 e 2.0; tronco vicino; superficie spostata di 4 e 10 cm. 3 seed.
#
# PREVISIONI SCRITTE PRIMA DELL'ESECUZIONE
#  P_Y1: nelle 8 condizioni ideali F ha |errore| mediano <= 1.5% e tutte le 24 esecuzioni danno una stima valida (>= 5 cm).
#  P_Y2: negli scenari con interni, guscio e superficie spostata, l'errore mediano di F differisce di al massimo 3 punti da quello
#        ottenuto prima delle correzioni (Tabella 1 del manoscritto: -5.7, -11.7, +11.5, +24.3, +0.2, +4.4): cioe' le correzioni
#        non cambiano cio' che mostrano i test su geometria ideale. Il tronco vicino e' esplorativo (prima: +24.5, vicino misurato
#        in 2 seed su 3).
#  P_Y3: il DBH diretto ha |errore| <= 3% in tutte le condizioni tranne "guscio 2.0", dove ci si aspetta circa +10%
#        (riferimento Python: ideali da -0.5 a 0.0; interni -1.3 e -2.6; guscio +2.8 e +9.6; vicino -0.3; spostata +1.5 e -0.2).
#
# NUOVA ESECUZIONE CON LA CONFIGURAZIONE FINALE (con la riparazione di connettivita' dei rami principali attiva).
# Le previsioni P_Y1-P_Y3 valevano per la prima esecuzione con il codice corretto (prima della riparazione; i valori sono in Tabella 1
# del manoscritto: F ideali da -0.1 a +0.4; interni 0.2 -5.4, interni 0.5 -10.7; guscio 0.5 +9.5, guscio 2.0 +6.9; vicino +24.9;
# spostata 4 cm +0.5, 10 cm +3.8; diretto invariato). PREVISIONI SCRITTE PRIMA DELLA NUOVA ESECUZIONE (08/10/2026):
#  P_Y4: la riparazione non cambia il DBH: in ognuna delle 15 condizioni l'errore mediano di F differisce di al massimo 1.0 punto dal
#        valore di Tabella 1 (nelle 8 ideali: tra -1.0 e +1.4).
#  P_Y5: nelle 8 condizioni ideali tutte le 24 esecuzioni danno una stima valida (>= 5 cm).
#  P_Y6: il DBH diretto e' deterministico: differisce da quello della prima esecuzione di al massimo 0.1 punti in ogni condizione.
# Prerequisiti (configurazione finale): surface_coverage.R corretto e ricaricato con source(), dbh_diretto.R caricato e
#   source(".../prova_tree_sets_fedele.R"); assign("define_main_branches", define_main_branches_fedele, envir = globalenv())
# Prerequisiti (sessione ripristinata): surface_coverage.R corretto e ricaricato con source(), dbh_diretto.R caricato.
# Uso:   source(".../prove_sintetiche.R");   r <- prove_sintetiche()          # circa 30-45 minuti
# ============================================================
prove_sintetiche <- function(cartella_R = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", file_out = NULL, offsets = 0:2) {
  mancanti <- c("treeqsm", "inputs_prod", "dbh_cylinder", "dbh_diretto", "surface_coverage"); mancanti <- mancanti[!sapply(mancanti, exists, envir = globalenv())]
  if (length(mancanti)) stop("mancano nella sessione: ", paste(mancanti, collapse = ", "), " (dbh_diretto: source('.../dbh_diretto.R'))", call. = FALSE)
  if (!exists("dbh_robusto", envir = globalenv())) cat("NOTA: dbh_robusto non e' nella sessione pulita: la colonna A (lettura personalizzata) sara' NA; non serve alla Tabella 1 finale\n")
  if (!("full" %in% names(formals(surface_coverage)))) stop("surface_coverage non e' la versione corretta: sostituisci il file e ricaricalo con source()")
  if (!exists("define_main_branches_fedele") || !identical(get("define_main_branches", envir = globalenv()), define_main_branches_fedele))
    stop("define_main_branches fedele non attiva: source('.../prova_tree_sets_fedele.R') e assign('define_main_branches', define_main_branches_fedele, envir = globalenv())")
  if (is.null(file_out)) file_out <- file.path(cartella_R, "prove_sintetiche_v5_riparazione.rds")
  
  rs_salvata <- get("relative_size", envir = globalenv())
  on.exit(assign("relative_size", rs_salvata, envir = globalenv()), add = TRUE)
  source(file.path(cartella_R, "relative_size.R"))
  if (!any(grepl("ChildSegment", deparse(get("relative_size", envir = globalenv()))))) stop("relative_size fedele non attiva")
  treeqsm_orig <- get("treeqsm", envir = globalenv())
  on.exit(assign("treeqsm", treeqsm_orig, envir = globalenv()), add = TRUE)
  txt <- readLines(file.path(cartella_R, "treeqsm.R"), warn = FALSE)
  i_td <- grep("treedata = tdata\\$treedata,", txt)
  if (length(i_td) != 1) stop("treeqsm.R: punto di inserimento non trovato")
  txt[i_td] <- paste0(txt[i_td], "\n            trunk_points = trunk_points,")
  tf <- tempfile(fileext = ".R"); writeLines(txt, tf); source(tf)
  assign(".off", 0L, envir = globalenv())
  assign("set.seed", function(seed, ...) base::set.seed(as.integer((as.numeric(seed) + get(".off", envir = globalenv()) * 7919) %% 2147483647), ...),
         envir = globalenv())
  on.exit(if (exists("set.seed", envir = globalenv(), inherits = FALSE)) rm("set.seed", envir = globalenv()), add = TRUE)
  
  gv <- function(x) if (is.null(x) || length(x) == 0 || !is.finite(x[1])) NA_real_ else as.numeric(x[1])
  dbh_fed <- function(m) {
    if (is.null(m$trunk_points)) return(NA_real_)
    cy <- m$cylinder
    r <- tryCatch(dbh_cylinder(cy$branch == 1, cy, seq_along(cy$radius), m$trunk_points), error = function(e) NULL)
    if (is.null(r)) NA_real_ else r$DBHcyl
  }
  dbh_catena_a <- function(m, zq) {
    cy <- m$cylinder; st <- which(cy$branch == 1)
    if (length(st) == 0) return(NA_real_)
    z0 <- cy$start[st, 3]; z1 <- z0 + cy$axis[st, 3] * cy$length[st]
    j <- which(pmin(z0, z1) <= zq & pmax(z0, z1) >= zq)
    if (length(j) == 0) return(NA_real_)
    200 * cy$radius[st[j[1]]]
  }
  gen <- function(r, L = 8, dens = 2000, sd = 0.005, cop = 1, f_int = 0, f_shell = 0, vicino = FALSE, shift = 0, seed = 1) {
    base::set.seed(seed)
    sup <- function(r, cx = 0, cy = 0) {
      n <- round(dens * 2 * pi * r * L)
      th <- runif(n, 0, 2 * pi * cop); z <- runif(n, 0, L); rr <- r + rnorm(n, 0, sd)
      cbind(cx + rr * cos(th), cy + rr * sin(th), z)
    }
    P <- sup(r); n_sup <- nrow(P)
    if (f_int > 0)   { n <- round(f_int * n_sup);   rr <- r * sqrt(runif(n)); th <- runif(n, 0, 2 * pi); P <- rbind(P, cbind(rr * cos(th), rr * sin(th), runif(n, 0, L))) }
    if (f_shell > 0) { n <- round(f_shell * n_sup); rr <- r + runif(n, 0, 0.15); th <- runif(n, 0, 2 * pi); P <- rbind(P, cbind(rr * cos(th), rr * sin(th), runif(n, 0, 3))) }
    if (vicino) P <- rbind(P, sup(0.15, cx = 0.9, cy = 0))
    if (shift > 0) { S <- sup(r); S[, 1] <- S[, 1] + shift; P <- rbind(P, S) }
    P
  }
  
  cond <- list()
  for (r in c(0.12, 0.25)) for (sdr in c(0.005, 0.010)) for (cp in c(1.0, 0.6))
    cond[[length(cond) + 1]] <- list(nome = sprintf("ideale r=%.2f rumore=%.3f copertura=%.1f", r, sdr, cp), tipo = "ideale", par = list(r = r, sd = sdr, cop = cp))
  cond <- c(cond, list(
    list(nome = "interni 0.2", tipo = "scenario", vecchio = -5.7, par = list(r = 0.12, f_int = 0.2)),
    list(nome = "interni 0.5", tipo = "scenario", vecchio = -11.7, par = list(r = 0.12, f_int = 0.5)),
    list(nome = "guscio 0.5", tipo = "scenario", vecchio = 11.5, par = list(r = 0.12, f_shell = 0.5)),
    list(nome = "guscio 2.0", tipo = "scenario", vecchio = 24.3, par = list(r = 0.12, f_shell = 2.0)),
    list(nome = "tronco vicino a 0.9 m", tipo = "vicino", vecchio = 24.5, par = list(r = 0.12, vicino = TRUE)),
    list(nome = "superficie spostata 4 cm", tipo = "scenario", vecchio = 0.2, par = list(r = 0.12, shift = 0.04)),
    list(nome = "superficie spostata 10 cm", tipo = "scenario", vecchio = 4.4, par = list(r = 0.12, shift = 0.10))))
  
  regs <- list()
  for (cd in cond) {
    for (o in offsets) {
      assign(".off", as.integer(o), envir = globalenv())
      P <- do.call(gen, c(cd$par, list(seed = 100 + o))); Pp <- P[P[, 3] >= 0.15, , drop = FALSE]
      q <- NULL; msg <- NA_character_
      invisible(capture.output(q <- tryCatch(treeqsm(Pp, inputs_prod), error = function(e) { msg <<- conditionMessage(e); NULL })))
      vero <- 200 * cd$par$r
      dd <- tryCatch(dbh_diretto(Pp)$diametro_cm, error = function(e) NA_real_)
      if (is.null(q) || length(q) == 0) {
        regs[[length(regs) + 1]] <- data.frame(cond = cd$nome, tipo = cd$tipo, vero = vero, off = o, A = NA, F = NA, E = NA, Diretto = dd, errore = msg, stringsAsFactors = FALSE)
        next
      }
      fcm <- 100 * vapply(q, function(m) gv(dbh_fed(m)), numeric(1)); ecm <- vapply(q, function(m) gv(dbh_catena_a(m, 1.3)), numeric(1))
      regs[[length(regs) + 1]] <- data.frame(cond = cd$nome, tipo = cd$tipo, vero = vero, off = o,
                                             A = if (exists("dbh_robusto", envir = globalenv())) tryCatch(dbh_robusto(q), error = function(e) NA_real_) else NA_real_,
                                             F = if (any(fcm >= 5, na.rm = TRUE)) median(fcm[fcm >= 5], na.rm = TRUE) else NA_real_,
                                             E = if (any(ecm >= 5, na.rm = TRUE)) median(ecm[ecm >= 5], na.rm = TRUE) else NA_real_,
                                             Diretto = dd, errore = NA_character_, stringsAsFactors = FALSE)
    }
    cat(sprintf("[%s] %s fatto\n", format(Sys.time(), "%H:%M:%S"), cd$nome))
    saveRDS(regs, file_out)
  }
  res <- do.call(rbind, regs)
  err <- function(x, v) 100 * (x - v) / v
  tb <- do.call(rbind, lapply(cond, function(cd) {
    g <- res[res$cond == cd$nome, ]; v <- g$vero[1]
    data.frame(condizione = cd$nome, vero_cm = v, errA = round(err(median(g$A, na.rm = TRUE), v), 1), errF = round(err(median(g$F, na.rm = TRUE), v), 1),
               errE = round(err(median(g$E, na.rm = TRUE), v), 1), errDiretto = round(err(median(g$Diretto, na.rm = TRUE), v), 1),
               F_prima_correzioni = if (is.null(cd$vecchio)) NA_real_ else cd$vecchio, stime_F_valide = sum(is.finite(g$F)), stringsAsFactors = FALSE)
  }))
  print(tb, row.names = FALSE)
  id <- tb[seq_len(8), ]; sc <- tb[9:15, ]
  cat(sprintf("\nP_Y1: nelle 8 condizioni ideali, massimo |errore| di F = %.2f%% (soglia 1.5%%); esecuzioni con stima F valida: %d su 24\n", max(abs(id$errF), na.rm = TRUE), sum(id$stime_F_valide)))
  cs <- sc[sc$condizione != "tronco vicino a 0.9 m", ]
  cat(sprintf("P_Y2: massimo scarto di F dai valori precedenti alle correzioni (6 scenari, escluso il tronco vicino) = %.1f punti (soglia 3)\n", max(abs(cs$errF - cs$F_prima_correzioni), na.rm = TRUE)))
  cat(sprintf("      tronco vicino: F = %+.1f%% (prima +24.5%%)\n", sc$errF[sc$condizione == "tronco vicino a 0.9 m"]))
  cat(sprintf("P_Y3: DBH diretto: massimo |errore| escluso 'guscio 2.0' = %.1f%% (soglia 3%%); guscio 2.0 = %+.1f%% (atteso circa +10%%)\n",
              max(abs(tb$errDiretto[tb$condizione != "guscio 2.0"]), na.rm = TRUE), tb$errDiretto[tb$condizione == "guscio 2.0"]))
  invisible(list(tabella = tb, esecuzioni = res))
}