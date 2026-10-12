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

# DIAGNOSI_STADI_R: stessi passi e stesse colonne di diagnosi_stadi_octave.m, per il porting in R (sessione pulita gia' caricata).
# Uso:  source(".../diagnosi_stadi_R.R");  d35 <- diagnosi_stadi_R(35)
diagnosi_stadi_R <- function(N = 35, cartella = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti", seed = 1) {
  stopifnot(exists("inputs_prod"), exists("cylinders"), exists("relative_size"))
  P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", N)))); mode(P) <- "numeric"; P <- P[P[, 3] >= 0.15, 1:3, drop = FALSE]
  inp <- inputs_prod; set.seed(seed); cat(sprintf("albero %d: %d punti\n", N, nrow(P)))
  out <- NULL
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
        cy <- cylinders(P, ts2$cover, s2, I)
        T <- which(cy$branch == 1); r <- cy$radius[T]; L <- cy$length[T]; cs <- cumsum(L); k <- which(cs >= 1.3)[1]
        d130 <- if (is.na(k)) NA_real_ else 200 * r[k]
        out <- rbind(out, data.frame(h = h, i = i, j = j, ns1 = length(ts1$cover$center), seg1 = length(s1$segments), lay1 = length(s1$segments[[1]]),
          RSmed = median(RS), RSmean = mean(RS), ns2 = length(ts2$cover$center), seg2 = length(s2$segments), lay2 = length(s2$segments[[1]]),
          ncil = length(cy$radius), ncilT = length(T), r3base = 100 * mean(r[seq_len(min(3, length(r)))]), D130cm = d130, Lfusto = sum(L)))
      }
    }
  }
  print(round(out, 2), row.names = FALSE)
  cat("\nMediane sugli 8 modelli:\n"); print(round(sapply(out[, -(1:3)], median, na.rm = TRUE), 2))
  invisible(out)
}
