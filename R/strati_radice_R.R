# STRATI_RADICE_R: struttura degli strati del segmento radice (il fusto) nel porting; stessi numeri di strati_radice_octave.m.
# Uso:  source(".../strati_radice_R.R");  s <- strati_radice_R(35)
strati_radice_R <- function(N = 35, cartella = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti", seed = 1) {
  P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", N)))); mode(P) <- "numeric"; P <- P[P[, 3] >= 0.15, 1:3, drop = FALSE]
  inp <- inputs_prod; set.seed(seed); zq <- c(0.25, 0.5, 1, 1.5, 2); res <- NULL; np10 <- NULL; ns10 <- NULL; info <- NULL
  for (h in seq_along(inp$PatchDiam1)) {
    I <- inp; I$PatchDiam1 <- inp$PatchDiam1[h]; I$BallRad1 <- inp$BallRad1[h]
    c1 <- cover_sets(P, I); ts1 <- tree_sets(P, c1, I); s1 <- segments(ts1$cover, ts1$Base, ts1$Forb); s1 <- correct_segments(P, ts1$cover, s1, I, 0, 1, 1)
    for (i in seq_along(inp$PatchDiam2Max)) { I$PatchDiam2Max <- inp$PatchDiam2Max[i]; I$BallRad2 <- inp$BallRad2[i]
      for (j in seq_along(inp$PatchDiam2Min)) { I$PatchDiam2Min <- inp$PatchDiam2Min[j]
        RS <- relative_size(P, ts1$cover, s1); c2 <- cover_sets(P, I, RS); ts2 <- tree_sets(P, c2, I, s1)
        Forb2 <- if (is.logical(ts2$Forb)) which(ts2$Forb) else ts2$Forb
        s2 <- segments(ts2$cover, ts2$Base, Forb2); s2 <- correct_segments(P, ts2$cover, s2, I, 1, 1, 0)
        SPar <- s2$ParentSegment; if (is.null(dim(SPar))) SPar <- cbind(as.integer(SPar), 0L); root <- which(SPar[, 1] == 0L)[1]
        Seg <- s2$segments[[root]]; nl <- length(Seg); cen <- as.integer(ts2$cover$center)
        nset <- sapply(Seg, length); zc <- sapply(Seg, function(S) mean(P[cen[unlist(S)], 3]))
        npt <- sapply(Seg, function(S) sum(sapply(unlist(S), function(s) length(ts2$cover$ball[[s]]))))
        res <- rbind(res, sapply(zq, function(z) { k <- which(zc >= z)[1]; if (is.na(k)) NA else k }))
        np10 <- rbind(np10, npt[1:10]); ns10 <- rbind(ns10, nset[1:10]); info <- rbind(info, c(nl, sum(nset), sum(npt), zc[1], max(zc)))
      } } }
  colnames(res) <- paste0("strato_a_z", zq); colnames(info) <- c("strati", "insiemi", "punti", "z_primo_strato", "z_ultimo")
  cat("\nStrato (indice) in cui la quota media dei centri raggiunge z, per ciascun modello:\n"); print(res, row.names = FALSE)
  cat("\nStrati, insiemi, punti del segmento radice; quota del primo e dell'ultimo strato:\n"); print(round(info, 2), row.names = FALSE)
  cat("\nMEDIANE: strati a z:", apply(res, 2, median, na.rm = TRUE), "| punti nei primi 10 strati:", apply(np10, 2, median), "| insiemi nei primi 10 strati:", apply(ns10, 2, median), "\n")
  invisible(list(res = res, info = info, np10 = np10, ns10 = ns10))
}
