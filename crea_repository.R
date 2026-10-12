# ============================================================
# CREA IL CONTENUTO DEL REPOSITORY a partire dalla cartella R validata (stato "v0.2.0-submission": configurazione FEDELE; file copiati come sono, senza modifiche)
# Prima di lanciarlo:
#  0. in una sessione R PULITA: source(".../sessione_pulita.R"); carica_libreria_pulita() (fedele = TRUE) per riscrivere manifesto_sessione_pulita.csv con la configurazione fedele
#  1. salvare una copia del manifesto della libreria validata (il file viene riscritto a ogni prova):
#       file.copy("F:/.../R/manifesto_sessione_pulita.csv", "F:/.../R/manifesto_validato.csv")
#  2. decomprimere lo scheletro (README.md, .gitignore, CITATION.cff, NOTICE.md, docs/, python/) in una cartella nuova, per esempio F:/EB/tartarino/TreeQSM-R
# Uso:  source("crea_repository.R");  crea_repository()
# Cosa fa: copia in <destinazione>/R i file della libreria (elenco esplicito del caricatore), create_input.R, la riparazione di connettivita' fedele e gli script dell'analisi;
# copia il manifesto in docs/manifesto_validato.csv; scrive docs/file_md5.txt; segnala i file mancanti, i percorsi assoluti e i file senza intestazione TREEQSM/GPL.
# NON copia dati, RDS, snapshot, .RData, ne' gli script operativi.
# ============================================================
crea_repository <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", destinazione = "F:/EB/tartarino/TreeQSM-R", includi_risultati_octave = FALSE, octave_dir = "D:/octave_treeqsm/TreeQSM-master/src") {
  libreria <- c("adjustments.R", "branches.R", "connected_components.R", "correct_segments.R", "cover_sets.R", "cross_product.R", "cubical_partition.R", "cylinder_fitting.R", "cylinders.R", "define_main_branches.R", "distances_between_lines.R", "distances_to_line.R", "growth_volume_correction.R", "helpers.R", "average.R", "least_squares_circle_centre.R", "least_squares_fitting.R", "make_tree_connected.R", "mat_vec_subtraction.R", "modify_parent.R", "modify_topology.R", "norm_vec.R", "parent_cylinder.R", "point_model_distance.R", "relative_size.R", "remove_small.R", "search_branch_top.R", "search_stem_top.R", "segment_direction.R", "segments.R", "select_cylinders.R", "surface_coverage.R", "surface_coverage2.R", "surface_coverage_filtering.R", "tree_data.R", "tree_data_core.R", "tree_data_distributions.R", "tree_sets.R", "treeqsm.R", "update_tree_data.R", "verticalcat.R", "correzioni_sessione_oggi.R", "create_input.R", "least_squares_cylinder_fedele.R")
  config <- c("prova_tree_sets_fedele.R")
  script <- c("raccolta_modelli_quote.R", "sessione_pulita.R", "lancia_pulita.R", "dbh_diretto.R", "valida_dbh_diretto.R", "dbh_diretto_quote.R", "esegui_test_r.R",
              "prove_sintetiche.R", "diagnosi_errori_silenziosi.R", "estrazione_alberi.R", "verifica_estrazione.R",
              "analisi_fedele.R", "confronta_octave_fedele.R", "tabella7_fedele.R", "riproduci_fedele.R", "prova_base_fedele.R", "prova_gn_fedele.R",
              "diagnosi_stadi_R.R", "profilo_fusto_R.R", "strati_radice_R.R", "fit_radice_R.R", "variabilita_base_casuale.R", "prova_orthonormal_casuale.R")
  if (!dir.exists(destinazione)) stop("La cartella di destinazione non esiste: decomprimere prima lo scheletro in ", destinazione)
  dest_R <- file.path(destinazione, "R"); dir.create(dest_R, showWarnings = FALSE); dir.create(file.path(destinazione, "docs"), showWarnings = FALSE)
  tutti <- c(libreria, config, script); presenti <- file.exists(file.path(base, tutti))
  if (any(!presenti)) cat("NON TROVATI in", base, ":", paste(tutti[!presenti], collapse = ", "), "\n")
  ok <- file.copy(file.path(base, tutti[presenti]), dest_R, overwrite = TRUE); cat(sprintf("copiati %d file su %d in %s\n", sum(ok), length(tutti), dest_R))
  man <- file.path(base, "manifesto_validato.csv")
  if (file.exists(man)) { file.copy(man, file.path(destinazione, "docs", "manifesto_validato.csv"), overwrite = TRUE); cat("manifesto validato copiato in docs/\n") }
  else cat("ATTENZIONE: manca manifesto_validato.csv in", base, "(vedere le istruzioni in testa a questo file)\n")
  si <- file.path(base, "docs_sessionInfo_fedele.txt")
  if (file.exists(si)) { file.copy(si, file.path(destinazione, "docs", "sessionInfo.txt"), overwrite = TRUE); cat("sessionInfo copiato in docs/sessionInfo.txt\n") } else cat("ATTENZIONE: manca docs_sessionInfo_fedele.txt in", base, "(si scrive con riproduci_fedele)\n")
  if (isTRUE(includi_risultati_octave)) {
    dir.create(file.path(destinazione, "results"), showWarnings = FALSE)
    for (cs in list.files(octave_dir, pattern = "^risultati_octave_.*\\.csv$", full.names = TRUE)) file.copy(cs, file.path(destinazione, "results", basename(cs)), overwrite = TRUE)
    cat("risultati di Octave copiati in results/ (decisione degli autori: sono stime per albero, non i dati grezzi)\n")
  }
  f <- list.files(dest_R, full.names = TRUE); f <- f[basename(f) != ".gitkeep"]
  writeLines(sprintf("%s  %s", unname(tools::md5sum(f)), basename(f)), file.path(destinazione, "docs", "file_md5.txt"))
  ass <- vapply(f, function(p) sum(grepl("[A-Za-z]:[/\\\\]", readLines(p, warn = FALSE))), numeric(1)); ass <- ass[ass > 0]
  cat(sprintf("\nfile con percorsi assoluti (da parametrizzare in un secondo momento): %d\n", length(ass)))
  if (length(ass) > 0) print(setNames(as.integer(ass), basename(names(ass))))
  senza <- basename(f)[!vapply(f, function(p) any(grepl("part of (the R port of )?TREEQSM", readLines(p, n = 8, warn = FALSE), ignore.case = TRUE)), logical(1))]
  cat(sprintf("\nfile senza intestazione 'This file is part of TREEQSM' nelle prime righe: %d\n", length(senza))); if (length(senza) > 0) cat(paste(senza, collapse = ", "), "\n")
  cat("\nModello di intestazione: docs/INTESTAZIONE_GPL.txt (da aggiungere in un secondo momento ai file derivati da TREEQSM; sono commenti e non cambiano il manifesto).\n")
  invisible(list(copiati = tutti[presenti], mancanti = tutti[!presenti], percorsi_assoluti = names(ass), senza_intestazione = senza))
}
