# PROVA SUI 21 ALBERI: base fedele del secondo passaggio (define_base_forb_fedele), con o senza Gauss-Newton, contro l'originale in Octave.
# Uso (sessione pulita con carica_libreria_pulita() eseguita; serve anche  source(".../prova_gn_fedele.R")  se gn = TRUE):
#   source(".../prova_base_fedele.R")
#   r <- prova_base_fedele("F:/.../risultati_octave_alberi_corretti_v2.csv")                 # solo base fedele; circa 12-15 minuti
#   r2 <- prova_base_fedele("F:/.../risultati_octave_alberi_corretti_v2.csv", gn = TRUE)     # base fedele + Gauss-Newton
# Non modifica nulla su disco: le funzioni sostituite vengono ripristinate alla fine.
# PREVISIONI SCRITTE PRIMA: P_H1 la mediana |errore| del DBH resta tra 2.8% e 4.2%; P_H2 la differenza mediana assoluta dall'originale (per albero)
#  scende da quella attuale (porting di produzione) di almeno 0.3 cm; P_H3 gli alberi 4, 16, 35 si avvicinano all'originale (|differenza| <= 1.5 cm in almeno 2 dei 3).
prova_base_fedele <- function(csv_octave, gn = FALSE, alberi = NULL, cartella_nuvola = "alberi_corretti", base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R") {
  stopifnot(exists("define_base_forb_fedele"), exists("raccolta_modelli_quote"), file.exists(csv_octave))
  if (gn) stopifnot(exists("least_squares_cylinder_gn"))
  ob <- get("define_base_forb", envir = globalenv()); ol <- get("least_squares_cylinder", envir = globalenv())
  on.exit({ assign("define_base_forb", ob, envir = globalenv()); assign("least_squares_cylinder", ol, envir = globalenv()) }, add = TRUE)
  assign("define_base_forb", define_base_forb_fedele, envir = globalenv()); if (gn) assign("least_squares_cylinder", least_squares_cylinder_gn, envir = globalenv())
  cart <- file.path(base, cartella_nuvola); tmp <- tempfile(fileext = ".rds"); t0 <- Sys.time()
  nu <- raccolta_modelli_quote(cart, cartella_R = base, file_out = tmp, offsets = 0, alberi = alberi)$modelli
  cat(sprintf("\nTempo: %.1f minuti\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rif <- read.csv(file.path(cart, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
  pr <- readRDS(file.path(cart, "modelli_quote_vol.rds"))$modelli; pr <- pr[pr$off == 0, ]; o <- read.csv(csv_octave)
  med <- function(x) { x <- x[is.finite(x) & x >= 5]; if (length(x)) median(x) else NA_real_ }
  ns <- sort(unique(nu$N))
  tb <- do.call(rbind, lapply(ns, function(n) { v <- rif$ref_cm[rif$N == n]
    data.frame(N = n, ref = v, Oct = med(o$dbh_cyl_cm[o$N == n]), prod = med(pr$dbh_matlab_cm[pr$N == n]), nuovo = med(nu$dbh_matlab_cm[nu$N == n])) }))
  tb$errOct <- 100 * (tb$Oct - tb$ref) / tb$ref; tb$errProd <- 100 * (tb$prod - tb$ref) / tb$ref; tb$errNuovo <- 100 * (tb$nuovo - tb$ref) / tb$ref
  tb$dProd <- tb$prod - tb$Oct; tb$dNuovo <- tb$nuovo - tb$Oct
  cat(sprintf("\n=== DBH (cm), mediana dei modelli, un seed; 'nuovo' = base fedele%s ===\n", if (gn) " + Gauss-Newton" else "")); print(round(tb, 1), row.names = FALSE)
  s <- function(x) c(med_abs = median(abs(x), na.rm = TRUE), media_abs = mean(abs(x), na.rm = TRUE), segno = median(x, na.rm = TRUE))
  cat("\nErrore %:  originale (Octave) | porting di produzione | porting con base fedele\n"); print(round(rbind(Octave = s(tb$errOct), produzione = s(tb$errProd), nuovo = s(tb$errNuovo)), 2))
  cat(sprintf("\nDifferenza assoluta dall'originale per albero (cm): mediana produzione %.2f | nuovo %.2f ; alberi entro 1 cm: produzione %d | nuovo %d su %d\n",
              median(abs(tb$dProd), na.rm = TRUE), median(abs(tb$dNuovo), na.rm = TRUE), sum(abs(tb$dProd) <= 1, na.rm = TRUE), sum(abs(tb$dNuovo) <= 1, na.rm = TRUE), nrow(tb)))
  cat(sprintf("Correlazione degli errori con l'originale: produzione %.2f | nuovo %.2f\n", cor(tb$errProd, tb$errOct, use = "complete.obs"), cor(tb$errNuovo, tb$errOct, use = "complete.obs")))
  invisible(list(tabella = tb, modelli = nu))
}
