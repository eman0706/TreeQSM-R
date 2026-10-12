# ============================================================
# PROVA DA SESSIONE PULITA
# Obiettivo: verificare che i risultati del lavoro si riproducano da una sessione R NUOVA, caricando la libreria dai file e applicando solo la
# procedura documentata, SENZA lo snapshot (snapshot_sessione.RData) e senza oggetti rimasti in memoria da sessioni precedenti.
# Procedura nel codice: caricamento ESPLICITO dei 41 file della libreria (piu' connected_components.R e create_input.R) che definiscono le funzioni raggiunte da treeqsm (elenco ricavato dal grafo delle chiamate di
# Figura 6; nessun bootstrap, nessuna esclusione da mantenere), poi surface_coverage.R e update_tree_data.R corretti (cartella correzioni_sorgenti), riparazione di
# connettivita' fedele (prova_tree_sets_fedele.R e assign), inputs_prod, raccolta_modelli_quote.R (che ri-carica da se' relative_size.R fedele e applica la patch a treeqsm).
# PERCHE' NON PIU' IL BOOTSTRAP: il primo tentativo con bootstrap_sessione.R ha eseguito tre script operativi non esclusi (param_test.R, param_tuning.R, test_nuovo.R: setwd,
# lettura del LAS, treeqsm, scrittura di risultati) e ha riscritto risultati/param_tuning.csv. Con l'elenco esplicito cambiano solo due definizioni in vigore rispetto al
# caricamento completo: relative_size (diventa la versione fedele, come nelle corse validate, che la ri-caricavano) e set.seed (resta quella di base: i due script di validazione
# la mascheravano con una versione che dipende dalla variabile globale .off).
# PREVISIONI SCRITTE PRIMA
#  P_C1: il caricamento esplicito termina senza errori (nessuna funzione "non trovata"): lo snapshot e gli altri 70 file non servono.
#  P_C2: in sessione pulita, per gli alberi 4, 16 e 35 (seed 0) su entrambe le nuvole, gli 8 modelli per albero riproducono tutte le grandezze registrate nelle corse
#        precedenti (DBHcyl ufficiale, diametri alle tre quote, volumi del fusto) entro 0.05 cm per i diametri e 0.1 litri per i volumi.
#  P_C3: gli oggetti che non sono funzioni nell'ambiente globale sono solo quelli creati dal caricatore.
#  P_C4: nessun file della cartella viene modificato oltre al manifesto (controllo con le date di modifica prima e dopo).
#  P_C5: set.seed e' quella di base e relative_size e' la versione fedele (usa ChildSegment).
#  P_C6: average e' la variante colMeans di average.R. Scoperta della prova: helpers.R e average.R definiscono average in due modi matematicamente equivalenti (colSums(X)/n e colMeans(P)) che
#        differiscono nell'ultimo bit; le corse validate usavano colMeans (versione dello snapshot): con colSums/n, un modello su 48 (albero 4, modello 5, nuvola 15) cambia il volume
#        del fusto sopra i 4 m di 1.35 L (0.2%); con colMeans tutto e' riprodotto esattamente. L'elenco di caricamento mette average.R dopo helpers.R perche' vinca colMeans.
# USO (R appena riavviato, senza aver caricato nulla; verificare che ls() sia vuoto):
#   source("F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/sessione_pulita.R")
#   carica_libreria_pulita()            # 1-2 minuti
#   r <- riproduci_pulita()             # circa 10 minuti (3 alberi, seed 0, due nuvole); riproduci_pulita(completa = TRUE) per tutti gli alberi (circa 45 minuti)
# Il manifesto manifesto_sessione_pulita.csv (funzione, file di definizione, impronta md5 del corpo) serve a confrontare in futuro il repository consolidato con questa libreria.
# ============================================================
carica_libreria_pulita <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", cartella_correzioni = NULL, fedele = TRUE) {
  altri <- setdiff(ls(globalenv(), all.names = TRUE), c("carica_libreria_pulita", "riproduci_pulita", "confronta_manifesto", ".Random.seed"))   # .Random.seed e' creato dall'uso del generatore, non e' un residuo
  if (length(altri) > 0) stop("La sessione non e' pulita: ", length(altri), " oggetti presenti (primi: ", paste(head(altri, 12), collapse = ", "), "). Probabile causa: R ha ripristinato un .RData all'avvio (cartella di lavoro: ", getwd(), "; .RData presente: ", file.exists(file.path(getwd(), ".RData")), "). Soluzione: RStudio > Tools > Global Options > General > togliere 'Restore .RData into workspace at startup' e riavviare R, oppure usare Rscript --vanilla. NON cancellare il .RData: rinominarlo.", call. = FALSE)
  options(keep.source = TRUE)
  prima <- file.info(list.files(base, recursive = TRUE, full.names = TRUE))
  if (!requireNamespace("lidR", quietly = TRUE)) cat("   nota: lidR non e' installato: serve solo per leggere i LAS e per l'estrazione degli alberi, non per riprodurre i modelli\n")
  elenco <- c("adjustments.R", "branches.R", "connected_components.R", "correct_segments.R", "cover_sets.R", "cross_product.R", "cubical_partition.R", "cylinder_fitting.R", "cylinders.R", "define_main_branches.R", "distances_between_lines.R", "distances_to_line.R", "growth_volume_correction.R", "helpers.R", "average.R", "least_squares_circle_centre.R", "least_squares_fitting.R", "make_tree_connected.R", "mat_vec_subtraction.R", "modify_parent.R", "modify_topology.R", "norm_vec.R", "parent_cylinder.R", "point_model_distance.R", "relative_size.R", "remove_small.R", "search_branch_top.R", "search_stem_top.R", "segment_direction.R", "segments.R", "select_cylinders.R", "surface_coverage.R", "surface_coverage2.R", "surface_coverage_filtering.R", "tree_data.R", "tree_data_core.R", "tree_data_distributions.R", "tree_sets.R", "treeqsm.R", "update_tree_data.R", "verticalcat.R", "correzioni_sessione_oggi.R")
  mancanti <- elenco[!file.exists(file.path(base, elenco))]; if (length(mancanti) > 0) stop("file della libreria non trovati in ", base, ": ", paste(mancanti, collapse = ", "))
  cat("1. caricamento esplicito di", length(elenco), "file della libreria\n")
  for (f in elenco) tryCatch(source(file.path(base, f)), error = function(e) stop("errore caricando ", f, ": ", conditionMessage(e), call. = FALSE))
  # create_input.R non e' raggiunto da treeqsm (i parametri si costruiscono fuori) ma serve per inputs_prod: contiene solo definizioni di funzione
  tryCatch(source(file.path(base, "create_input.R")), error = function(e) stop("errore caricando create_input.R: ", conditionMessage(e), call. = FALSE))
  cat("2. correzioni dei sorgenti (surface_coverage.R e update_tree_data.R)\n")
  cand <- c(cartella_correzioni, file.path(base, "correzioni_sorgenti"), file.path(dirname(base), "correzioni_sorgenti"), file.path(dirname(dirname(base)), "correzioni_sorgenti"))
  cand <- cand[!is.null(cand) & dir.exists(cand)]
  for (f in c("surface_coverage.R", "update_tree_data.R")) {
    p <- file.path(cand, f); p <- p[file.exists(p)]
    if (length(p) > 0) { source(p[1]); cat("   caricato", p[1], "\n") } else cat("   ATTENZIONE: non trovo", f, "nelle cartelle correzioni_sorgenti: si usa la versione dei moduli\n")
  }
  if (exists("crown_measures", envir = globalenv())) cat("   crown_measures: accesso a $U ancora presente (versione non corretta):", any(grepl("\\$U\\b", deparse(get("crown_measures", envir = globalenv())))), "(non e' sul percorso del DBH)\n")
  if (!("full" %in% names(formals(get("surface_coverage", envir = globalenv()))))) stop("surface_coverage non e' la versione corretta (manca l'argomento full): copiare correzioni_sorgenti/surface_coverage.R e passare cartella_correzioni")
  cat("3. riparazione di connettivita' fedele\n")
  source(file.path(base, "prova_tree_sets_fedele.R")); assign("define_main_branches", define_main_branches_fedele, envir = globalenv())
  stopifnot(identical(get("define_main_branches", envir = globalenv()), define_main_branches_fedele))
  if (fedele) {   # versione fedele a TreeQSM 2.4.1: base del secondo passaggio e fit del cilindro (Gauss-Newton) come nell'originale
    if (!file.exists(file.path(base, "least_squares_cylinder_fedele.R"))) stop("manca least_squares_cylinder_fedele.R in ", base)
    source(file.path(base, "least_squares_cylinder_fedele.R"))
    assign("define_base_forb", define_base_forb_fedele, envir = globalenv()); assign("least_squares_cylinder", least_squares_cylinder_fedele, envir = globalenv())
    stopifnot(identical(get("define_base_forb", envir = globalenv()), define_base_forb_fedele), identical(get("least_squares_cylinder", envir = globalenv()), least_squares_cylinder_fedele))
    cat("   attive: define_base_forb_fedele e least_squares_cylinder_fedele (Gauss-Newton)\n")
  } else cat("   ATTENZIONE: versione di produzione precedente (base a 5 insiemi, Levenberg-Marquardt)\n")
  cat("4. parametri di produzione\n")
  ip <- create_input(); ip$PatchDiam1 <- c(0.08, 0.12); ip$PatchDiam2Min <- c(0.02, 0.03); ip$PatchDiam2Max <- c(0.07, 0.1)
  ip$BallRad1 <- ip$PatchDiam1 + 0.015; ip$BallRad2 <- ip$PatchDiam2Max + 0.01; ip$Dist <- TRUE; assign("inputs_prod", ip, envir = globalenv())
  cat("5. raccolta_modelli_quote.R\n"); source(file.path(base, "raccolta_modelli_quote.R"))
  stopifnot(exists("treeqsm"), exists("inputs_prod"), exists("raccolta_modelli_quote"), exists("dbh_cylinder"))
  # P_C3: oggetti che non sono funzioni
  nomi <- ls(globalenv(), all.names = TRUE); nf <- nomi[!vapply(nomi, function(n) is.function(get(n, envir = globalenv())), logical(1))]
  cat("\nOggetti che non sono funzioni nell'ambiente globale (", length(nf), "): ", paste(nf, collapse = ", "), "\n", sep = "")
  # manifesto delle funzioni raggiunte da treeqsm
  visti <- character(0); coda <- "treeqsm"
  while (length(coda) > 0) {
    x <- coda[1]; coda <- coda[-1]
    if (x %in% visti || !exists(x, envir = globalenv(), inherits = FALSE)) next
    fx <- get(x, envir = globalenv()); if (!is.function(fx) || is.primitive(fx)) next
    visti <- c(visti, x); coda <- c(coda, intersect(all.names(body(fx)), nomi))
  }
  md5 <- function(fx) { f <- tempfile(); writeLines(deparse(fx), f); unname(tools::md5sum(f)) }
  man <- data.frame(funzione = visti, file = vapply(visti, function(x) { s <- utils::getSrcFilename(get(x, envir = globalenv())); if (length(s) == 0) NA_character_ else s[1] }, character(1)),
                    md5 = vapply(visti, function(x) md5(get(x, envir = globalenv())), character(1)), stringsAsFactors = FALSE)
  write.csv(man, file.path(base, "manifesto_sessione_pulita.csv"), row.names = FALSE)
  cat(sprintf("Manifesto: %d funzioni raggiunte da treeqsm (%d senza file di definizione noto) scritto in manifesto_sessione_pulita.csv\n", nrow(man), sum(is.na(man$file))))
  cat("P_C1: caricamento completato senza errori.\n")
  cat("P_C6: average e' la variante colMeans (quella delle corse validate):", any(grepl("colMeans", deparse(get("average", envir = globalenv())))), "\n")
  cat("P_C5: set.seed di base:", !exists("set.seed", envir = globalenv(), inherits = FALSE), "| relative_size fedele (ChildSegment):", any(grepl("ChildSegment", deparse(get("relative_size", envir = globalenv())))), "\n")
  dopo <- file.info(list.files(base, recursive = TRUE, full.names = TRUE)); nuovi <- setdiff(rownames(dopo), rownames(prima))
  com <- intersect(rownames(dopo), rownames(prima)); mod <- com[dopo[com, "mtime"] != prima[com, "mtime"]]
  cat("P_C4: file nuovi o modificati dal caricamento:", if (length(c(nuovi, mod)) == 0) "nessuno" else paste(basename(c(nuovi, mod)), collapse = ", "), "(atteso: solo manifesto_sessione_pulita.csv)\n")
  invisible(man)
}
riproduci_pulita <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", alberi = c(4, 16, 35), completa = FALSE, nuvole = c("alberi_corretti", "alberi_corretti_las11")) {
  stopifnot(exists("raccolta_modelli_quote"), exists("define_main_branches_fedele"))
  col_d <- c("dbh_matlab_cm", "q130_qsm", "q130_cyl", "q200_qsm", "q200_cyl", "q400_qsm", "q400_cyl"); col_v <- c("v_13_40", "v_0_4", "v_0_8", "v_stem")
  esiti <- list(); confronti <- list()
  for (cl in nuvole) {
    cart <- file.path(base, cl); tmp <- tempfile(fileext = ".rds")
    nuovo <- raccolta_modelli_quote(cart, cartella_R = base, file_out = tmp, offsets = 0, alberi = if (completa) NULL else alberi)
    if (is.null(nuovo$modelli) || nrow(nuovo$modelli) == 0) stop("raccolta_modelli_quote non ha prodotto nessun modello su ", cl, ": treeqsm ha dato errore (vedere i messaggi ATTENZIONE sopra)", call. = FALSE)
    rif <- readRDS(file.path(cart, "modelli_quote_vol.rds"))$modelli; nu <- nuovo$modelli
    mm <- merge(nu, rif, by = c("N", "off", "modello"), suffixes = c("_n", "_r")); confronti[[cl]] <- mm
    cat(sprintf("\n=== %s: %d modelli nuovi, %d confrontati con le corse precedenti ===\n", cl, nrow(nu), nrow(mm)))
    ok_tutto <- nrow(mm) == nrow(nu) && nrow(nu) > 0
    for (cc in c(col_d, col_v)) {
      a <- mm[[paste0(cc, "_n")]]; b <- mm[[paste0(cc, "_r")]]; d <- abs(a - b); tol <- if (cc %in% col_v) 0.1 else 0.05
      ent <- d <= tol | (is.na(a) & is.na(b)); ok <- mean(ent, na.rm = TRUE)
      cat(sprintf("  %-14s entro %.2f: %5.1f%% | differenza massima %.4f\n", cc, tol, 100 * ok, suppressWarnings(max(d, na.rm = TRUE)))); ok_tutto <- ok_tutto && ok == 1
      fuori <- which(!ent)
      if (length(fuori) > 0) { cat("    modelli fuori tolleranza:\n"); print(data.frame(albero = mm$N[fuori], modello = mm$modello[fuori], nuovo = round(a[fuori], 4), prima = round(b[fuori], 4), differenza = round(d[fuori], 4)), row.names = FALSE) }
    }
    esiti[[cl]] <- ok_tutto
  }
  cat(sprintf("\nP_C2 (tutte le grandezze riprodotte, entrambe le nuvole): %s\n", if (all(unlist(esiti))) "CONFERMATA" else "NON CONFERMATA"))
  invisible(list(esiti = esiti, confronto = confronti))
}

confronta_manifesto <- function(nuovo, riferimento) {
  a <- read.csv(nuovo, stringsAsFactors = FALSE); b <- read.csv(riferimento, stringsAsFactors = FALSE)
  m <- merge(a, b, by = "funzione", all = TRUE, suffixes = c("_nuovo", "_rif"))
  diversi <- m[is.na(m$md5_nuovo) | is.na(m$md5_rif) | m$md5_nuovo != m$md5_rif, ]
  cat(sprintf("funzioni: %d nel nuovo manifesto, %d nel riferimento; con md5 diverso o presenti in uno solo: %d\n", nrow(a), nrow(b), nrow(diversi)))
  if (nrow(diversi) > 0) print(diversi[, c("funzione", "file_nuovo", "file_rif")], row.names = FALSE)
  invisible(diversi)
}