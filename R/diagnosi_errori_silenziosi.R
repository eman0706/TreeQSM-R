# ============================================================
# DIAGNOSI: errori che il porting inghiotte senza fermarsi durante treeqsm
#
# La lettura riga per riga ha mostrato che diverse funzioni del porting catturano gli errori e proseguono, e che i nostri script di
# raccolta nascondevano i messaggi stampati con capture.output():
#  - cylinders.R: [CYL-ERR] (cylinder_fitting), [ADJ-ERR] (adjustments), [CYL-LOOP-ERR] (intero segmento), e
#    tryCatch(parent_cylinder(...), error = function(e) NULL) del tutto silenzioso: se fallisce, il cilindro resta senza genitore;
#  - cylinder_fitting.R: [CF-WARN] "regione fallita": la regione che da' errore INTERROMPE il fit del segmento (il fusto si ferma li');
#  - adjustments (effettiva): tryCatch(least_squares_circle_centre(...), error = function(e) NULL);
#  - modify_topology, search_stem_top, search_branch_top: warning() con ripieghi se il ciclo non converge o i rapporti non sono definiti.
# Lo script esegue treeqsm (1 seed) su alcuni alberi SENZA scartare i messaggi, conta le occorrenze e registra gli errori di
# parent_cylinder e di least_squares_circle_centre senza cambiarne il comportamento.
#
# PREVISIONI SCRITTE PRIMA
#  P_E1: nessuna regione fallita ne' segmento interrotto: [CF-WARN], [CYL-ERR], [CYL-LOOP-ERR] e [ADJ-ERR] valgono 0 in tutti gli alberi.
#  P_E2: gli errori di parent_cylinder sono meno del 5% delle chiamate; quelli di least_squares_circle_centre meno del 5%.
#  P_E3: nessun avviso (warning) da modify_topology, search_stem_top e search_branch_top.
#
# Prerequisito: stessa configurazione della validazione finale: surface_coverage corretta e define_main_branches fedele attivo, cioe'
#   source(".../prova_tree_sets_fedele.R"); assign("define_main_branches", define_main_branches_fedele, envir = globalenv())
# Uso:  source(".../diagnosi_errori_silenziosi.R");  d <- diagnosi_errori_silenziosi(".../alberi_corretti")      # circa 10 minuti
# ============================================================
diagnosi_errori_silenziosi <- function(cartella, alberi = c(4, 35, 16), offset = 0, cartella_R = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R") {
  stopifnot(exists("treeqsm"), exists("inputs_prod"), exists("parent_cylinder"), exists("least_squares_circle_centre"), exists("surface_coverage"))
  if (!("full" %in% names(formals(surface_coverage)))) stop("surface_coverage non e' la versione corretta")
  if (exists("define_main_branches_fedele") && !identical(get("define_main_branches", envir = globalenv()), define_main_branches_fedele))
    cat("ATTENZIONE: define_main_branches non e' la versione fedele (configurazione diversa dalla validazione finale)\n")
  orig <- list(parent_cylinder = get("parent_cylinder", envir = globalenv()), lsc = get("least_squares_circle_centre", envir = globalenv()),
               relative_size = get("relative_size", envir = globalenv()))
  on.exit({ assign("parent_cylinder", orig$parent_cylinder, envir = globalenv()); assign("least_squares_circle_centre", orig$lsc, envir = globalenv())
            assign("relative_size", orig$relative_size, envir = globalenv()) }, add = TRUE)
  source(file.path(cartella_R, "relative_size.R"))
  assign(".off", as.integer(offset), envir = globalenv())
  assign("set.seed", function(seed, ...) base::set.seed(as.integer((as.numeric(seed) + get(".off", envir = globalenv()) * 7919) %% 2147483647), ...), envir = globalenv())
  on.exit(if (exists("set.seed", envir = globalenv(), inherits = FALSE)) rm("set.seed", envir = globalenv()), add = TRUE)
  reg <- new.env()
  avvolgi <- function(nome, f) {
    reg[[nome]] <- list(n = 0L, err = 0L, msg = character(0))
    function(...) { reg[[nome]]$n <- reg[[nome]]$n + 1L
      tryCatch(f(...), error = function(e) { reg[[nome]]$err <- reg[[nome]]$err + 1L; reg[[nome]]$msg <- c(reg[[nome]]$msg, conditionMessage(e)); stop(e) }) }
  }
  assign("parent_cylinder", avvolgi("parent_cylinder", orig$parent_cylinder), envir = globalenv())
  assign("least_squares_circle_centre", avvolgi("least_squares_circle_centre", orig$lsc), envir = globalenv())
  esiti <- list()
  for (n in alberi) {
    reg$parent_cylinder <- list(n = 0L, err = 0L, msg = character(0)); reg$least_squares_circle_centre <- list(n = 0L, err = 0L, msg = character(0))
    P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", n)))); mode(P) <- "numeric"; Pp <- P[P[, 3] >= 0.15, , drop = FALSE]
    avvisi <- character(0); q <- NULL
    out <- capture.output(q <- withCallingHandlers(tryCatch(treeqsm(Pp, inputs_prod), error = function(e) { cat("ERRORE-TREEQSM:", conditionMessage(e), "\n"); NULL }),
                                                   warning = function(w) { avvisi <<- c(avvisi, conditionMessage(w)); invokeRestart("muffleWarning") }))
    cnt <- function(p) sum(grepl(p, out, fixed = TRUE))
    esiti[[length(esiti) + 1]] <- data.frame(N = n, modelli = length(q), CF_WARN = cnt("[CF-WARN]"), CYL_ERR = cnt("[CYL-ERR]"), CYL_LOOP_ERR = cnt("[CYL-LOOP-ERR]"), ADJ_ERR = cnt("[ADJ-ERR]"),
      parentcyl_chiamate = reg$parent_cylinder$n, parentcyl_errori = reg$parent_cylinder$err, lsc_chiamate = reg$least_squares_circle_centre$n, lsc_errori = reg$least_squares_circle_centre$err,
      avvisi = length(avvisi), avvisi_modify_topology_search = sum(grepl("modify_topology|search_stem_top|search_branch_top", avvisi)))
    cat(sprintf("albero %d fatto\n", n))
    if (reg$parent_cylinder$err > 0) { cat("  messaggi di errore di parent_cylinder:\n"); print(head(sort(table(reg$parent_cylinder$msg), decreasing = TRUE), 4)) }
    if (reg$least_squares_circle_centre$err > 0) { cat("  messaggi di errore di least_squares_circle_centre:\n"); print(head(sort(table(reg$least_squares_circle_centre$msg), decreasing = TRUE), 4)) }
    if (length(avvisi) > 0) { cat("  avvisi:\n"); print(head(sort(table(avvisi), decreasing = TRUE), 4)) }
    righe_err <- unique(grep("ERR|WARN", out, value = TRUE)); if (length(righe_err) > 0) { cat("  righe stampate con ERR o WARN (prime 6 diverse):\n"); print(head(righe_err, 6)) }
  }
  tb <- do.call(rbind, esiti); print(tb, row.names = FALSE)
  cat(sprintf("\nP_E1: regioni fallite + errori di cylinders/adjustments in totale: %d (atteso 0)\n", sum(tb$CF_WARN + tb$CYL_ERR + tb$CYL_LOOP_ERR + tb$ADJ_ERR)))
  cat(sprintf("P_E2: errori di parent_cylinder %.1f%% delle chiamate; di least_squares_circle_centre %.1f%% (soglia 5%%)\n",
              100 * sum(tb$parentcyl_errori) / max(1, sum(tb$parentcyl_chiamate)), 100 * sum(tb$lsc_errori) / max(1, sum(tb$lsc_chiamate))))
  cat(sprintf("P_E3: avvisi da modify_topology o search_*: %d (atteso 0)\n", sum(tb$avvisi_modify_topology_search)))
  invisible(tb)
}
