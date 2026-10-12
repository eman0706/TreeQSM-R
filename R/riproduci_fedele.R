# RIPRODUCIBILITA' DELLA VERSIONE FEDELE DA SESSIONE PULITA
# In una sessione R appena avviata (senza .RData ripristinato):
#   source(".../sessione_pulita.R"); carica_libreria_pulita(); source(".../riproduci_fedele.R"); r <- riproduci_fedele()
# Rieseguire treeqsm (alberi 4, 16, 35; 3 seed; entrambe le nuvole) e confrontare ogni grandezza con modelli_quote_fedele.rds (stesso codice, stessi seed:
# attesa riproduzione esatta come per la versione precedente, 48/48 modelli). Scrive anche docs_sessionInfo_fedele.txt nella cartella R con sessionInfo() e i tempi.
# PREVISIONE SCRITTA PRIMA: P_X1 tutte le grandezze riprodotte entro le tolleranze (0.05 cm per i diametri, 0.1 L per i volumi) su 48/48 modelli.
riproduci_fedele <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", alberi = c(4, 16, 35), nuvole = c("alberi_corretti", "alberi_corretti_las11")) {
  stopifnot(exists("raccolta_modelli_quote"), exists("define_main_branches_fedele"), exists("define_base_forb_fedele"), exists("least_squares_cylinder_fedele"))
  stopifnot(identical(get("define_base_forb", envir = globalenv()), define_base_forb_fedele), identical(get("least_squares_cylinder", envir = globalenv()), least_squares_cylinder_fedele))
  col_d <- c("dbh_matlab_cm", "q130_qsm", "q130_cyl", "q200_qsm", "q200_cyl", "q400_qsm", "q400_cyl"); col_v <- c("v_13_40", "v_0_4", "v_0_8", "v_stem")
  esiti <- list(); tempi <- numeric(0)
  for (cl in nuvole) {
    cart <- file.path(base, cl); tmp <- tempfile(fileext = ".rds"); t0 <- Sys.time()
    nuovo <- raccolta_modelli_quote(cart, cartella_R = base, file_out = tmp, offsets = 0:2, alberi = alberi)$modelli
    tempi[cl] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (is.null(nuovo) || nrow(nuovo) == 0) stop("nessun modello prodotto su ", cl, call. = FALSE)
    rif <- readRDS(file.path(cart, "modelli_quote_fedele.rds"))$modelli
    mm <- merge(nuovo, rif, by = c("N", "off", "modello"), suffixes = c("_n", "_r"))
    cat(sprintf("\n=== %s: %d modelli nuovi, %d confrontati con la corsa fedele salvata (%.1f minuti) ===\n", cl, nrow(nuovo), nrow(mm), tempi[cl]))
    ok_tutto <- nrow(mm) == nrow(nuovo) && nrow(nuovo) > 0
    for (cc in c(col_d, col_v)) {
      a <- mm[[paste0(cc, "_n")]]; b <- mm[[paste0(cc, "_r")]]; d <- abs(a - b); tol <- if (cc %in% col_v) 0.1 else 0.05
      ent <- d <= tol | (is.na(a) & is.na(b)); ok <- mean(ent, na.rm = TRUE)
      cat(sprintf("  %-14s entro %.2f: %5.1f%% | differenza massima %.4f\n", cc, tol, 100 * ok, suppressWarnings(max(d, na.rm = TRUE)))); ok_tutto <- ok_tutto && ok == 1
      fuori <- which(!ent); if (length(fuori) > 0) print(data.frame(albero = mm$N[fuori], seed = mm$off[fuori], modello = mm$modello[fuori], nuovo = round(a[fuori], 4), salvato = round(b[fuori], 4)), row.names = FALSE)
    }
    esiti[[cl]] <- ok_tutto
  }
  cat(sprintf("\nP_X1 (tutte le grandezze riprodotte, entrambe le nuvole, %d alberi x 3 seed): %s\n", length(alberi), if (all(unlist(esiti))) "CONFERMATA" else "NON CONFERMATA"))
  writeLines(c(sprintf("Corsa di riproducibilita' della versione fedele, %s", format(Sys.time())), sprintf("Tempi (minuti, %d alberi x 3 seed x 8 modelli): %s", length(alberi), paste(names(tempi), round(tempi, 1), collapse = "; ")), "", capture.output(print(sessionInfo()))), file.path(base, "docs_sessionInfo_fedele.txt"))
  cat("sessionInfo e tempi scritti in", file.path(base, "docs_sessionInfo_fedele.txt"), "\n")
  invisible(esiti)
}
