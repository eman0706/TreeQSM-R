# PROFILO_FUSTO_R: come profilo_fusto_octave.m, per il porting (sessione pulita gia' caricata).
# Uso:  source(".../profilo_fusto_R.R");  p35 <- profilo_fusto_R(35)
profilo_fusto_R <- function(N = 35, cartella = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti", seed = 1) {
  stopifnot(exists("inputs_prod"), exists("cylinders"), exists("relative_size"))
  P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", N)))); mode(P) <- "numeric"; P <- P[P[, 3] >= 0.15, 1:3, drop = FALSE]
  inp <- inputs_prod; set.seed(seed); quote <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 6)
  fin <- unm <- sc <- NULL; ncl <- NULL
  for (h in seq_along(inp$PatchDiam1)) {
    I <- inp; I$PatchDiam1 <- inp$PatchDiam1[h]; I$BallRad1 <- inp$BallRad1[h]
    c1 <- cover_sets(P, I); ts1 <- tree_sets(P, c1, I)
    s1 <- segments(ts1$cover, ts1$Base, ts1$Forb); s1 <- correct_segments(P, ts1$cover, s1, I, 0, 1, 1)
    for (i in seq_along(inp$PatchDiam2Max)) {
      I$PatchDiam2Max <- inp$PatchDiam2Max[i]; I$BallRad2 <- inp$BallRad2[i]
      for (j in seq_along(inp$PatchDiam2Min)) {
        I$PatchDiam2Min <- inp$PatchDiam2Min[j]
        RS <- relative_size(P, ts1$cover, s1); c2 <- cover_sets(P, I, RS); ts2 <- tree_sets(P, c2, I, s1)
        Forb2 <- if (is.logical(ts2$Forb)) which(ts2$Forb) else ts2$Forb
        s2 <- segments(ts2$cover, ts2$Base, Forb2); s2 <- correct_segments(P, ts2$cover, s2, I, 1, 1, 0)
        invisible(capture.output(cy <- cylinders(P, ts2$cover, s2, I)))
        T <- which(cy$branch == 1); L <- cy$length[T]; cs <- cumsum(L); r <- cy$radius[T]; u <- cy$UnModRadius[T]; s <- cy$SurfCov[T]
        f <- g <- w <- rep(NA_real_, length(quote))
        for (q in seq_along(quote)) { k <- which(cs >= quote[q])[1]; if (!is.na(k)) { f[q] <- 200 * r[k]; g[q] <- 200 * u[k]; w[q] <- s[k] } }
        fin <- rbind(fin, f); unm <- rbind(unm, g); sc <- rbind(sc, w); ncl <- rbind(ncl, c(length(T), sum(abs(r - u) > 1e-6), mean(L)))
      }
    }
  }
  colnames(fin) <- colnames(unm) <- colnames(sc) <- paste0("L", quote); colnames(ncl) <- c("cil_fusto", "modificati", "lung_media")
  cat("\nDiametro DOPO adjustments (cm), colonne = lunghezza catena (m)\n"); print(round(fin, 1), row.names = FALSE)
  cat("Diametro PRIMA di adjustments (UnModRadius, cm)\n"); print(round(unm, 1), row.names = FALSE)
  cat("SurfCov del cilindro\n"); print(round(sc, 2), row.names = FALSE)
  cat("Cilindri del fusto | modificati da adjustments | lunghezza media (m)\n"); print(round(ncl, 2), row.names = FALSE)
  cat("\nMEDIANE dopo:\n"); print(round(apply(fin, 2, median, na.rm = TRUE), 1)); cat("MEDIANE prima:\n"); print(round(apply(unm, 2, median, na.rm = TRUE), 1)); cat("MEDIANE SurfCov:\n"); print(round(apply(sc, 2, median, na.rm = TRUE), 2))
  invisible(list(dopo = fin, prima = unm, sc = sc, ncl = ncl))
}
