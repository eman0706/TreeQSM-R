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

# PROVA DI UNITA' di cylinder_fitting: stessi dati d'ingresso in R e in Octave, nessuna casualita' a monte.
# Passo 1 (R, sessione pulita gia' caricata):  source(".../fit_radice_R.R");  esporta_fit_radice(35)
#   esegue gli stadi fino al secondo segments per UN modello, prende il segmento radice (il fusto) e salva in cartella_out:
#   fit_P.csv, fit_Points.csv, fit_Ind.csv (gli argomenti di cylinder_fitting) e fit_R_out.csv (il risultato del porting).
# Passo 2 (Octave, stessa cartella): confronta_fit_radice_octave   -> fit_Oct_out.csv
# Passo 3 (R):  confronta_fit_radice()   -> tabella cilindro per cilindro dei due risultati.
esporta_fit_radice <- function(N = 35, h = 2, i = 1, j = 1, cartella = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/alberi_corretti", cartella_out = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/prova_fit", seed = 1) {
  dir.create(cartella_out, showWarnings = FALSE)
  P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", N)))); mode(P) <- "numeric"; P <- P[P[, 3] >= 0.15, 1:3, drop = FALSE]
  inp <- inputs_prod; set.seed(seed)
  I <- inp; I$PatchDiam1 <- inp$PatchDiam1[h]; I$BallRad1 <- inp$BallRad1[h]
  c1 <- cover_sets(P, I); ts1 <- tree_sets(P, c1, I); s1 <- segments(ts1$cover, ts1$Base, ts1$Forb); s1 <- correct_segments(P, ts1$cover, s1, I, 0, 1, 1)
  I$PatchDiam2Max <- inp$PatchDiam2Max[i]; I$BallRad2 <- inp$BallRad2[i]; I$PatchDiam2Min <- inp$PatchDiam2Min[j]
  RS <- relative_size(P, ts1$cover, s1); c2 <- cover_sets(P, I, RS); ts2 <- tree_sets(P, c2, I, s1)
  Forb2 <- if (is.logical(ts2$Forb)) which(ts2$Forb) else ts2$Forb
  s2 <- segments(ts2$cover, ts2$Base, Forb2); s2 <- correct_segments(P, ts2$cover, s2, I, 1, 1, 0)
  SPar <- s2$ParentSegment; if (is.null(dim(SPar))) SPar <- cbind(as.integer(SPar), 0L)
  root <- which(SPar[, 1] == 0L)[1]; Seg <- s2$segments[[root]]; nl <- length(Seg)
  vc <- verticalcat(Seg); Sets <- vc$Sets; IndSets <- vc$IndSets
  BallSize <- sapply(Sets, function(s) length(ts2$cover$ball[[s]])); Points <- unlist(lapply(Sets, function(s) ts2$cover$ball[[s]]))
  Ind <- matrix(1L, nl, 2); for (k in seq_len(nl)) Ind[k, 2] <- sum(BallSize[IndSets[k, 1]:IndSets[k, 2]])
  Ind[, 2] <- cumsum(Ind[, 2]); if (nl > 1) Ind[2:nl, 1] <- Ind[2:nl, 1] + Ind[1:(nl - 1), 2]
  cat(sprintf("segmento radice %d: %d strati, %d punti (con ripetizioni)\n", root, nl, length(Points)))
  write.csv(P, file.path(cartella_out, "fit_P.csv"), row.names = FALSE); write.csv(data.frame(Points = Points), file.path(cartella_out, "fit_Points.csv"), row.names = FALSE)
  write.csv(data.frame(a = Ind[, 1], b = Ind[, 2]), file.path(cartella_out, "fit_Ind.csv"), row.names = FALSE)
  fit <- cylinder_fitting(P, Points, Ind, nl, root); cy <- fit$cyl
  out <- data.frame(radius = cy$radius, length = cy$length, sx = cy$start[, 1], sy = cy$start[, 2], sz = cy$start[, 3], ax = cy$axis[, 1], ay = cy$axis[, 2], az = cy$axis[, 3],
                    SurfCov = cy$SurfCov, mad = cy$mad, conv = as.numeric(cy$conv), rel = as.numeric(cy$rel))
  write.csv(out, file.path(cartella_out, "fit_R_out.csv"), row.names = FALSE)
  cat(sprintf("salvati i file in %s: %d cilindri nel risultato del porting\n", cartella_out, nrow(out)))
  invisible(out)
}
confronta_fit_radice <- function(cartella_out = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R/prova_fit") {
  r <- read.csv(file.path(cartella_out, "fit_R_out.csv")); o <- read.csv(file.path(cartella_out, "fit_Oct_out.csv"), header = FALSE)
  names(o) <- c("radius", "length", "sx", "sy", "sz", "ax", "ay", "az", "SurfCov", "mad", "conv", "rel")
  cat(sprintf("cilindri: porting %d, originale %d\n", nrow(r), nrow(o)))
  f <- function(d) data.frame(z = round(d$sz, 3), L = round(d$length, 3), D_cm = round(200 * d$radius, 1), SC = round(d$SurfCov, 2), conv = d$conv, rel = d$rel)
  cat("\n=== PORTING (R) ===\n"); print(head(f(r), 25), row.names = FALSE); cat("\n=== ORIGINALE (Octave) ===\n"); print(head(f(o), 25), row.names = FALSE)
  qf <- function(d, q) { cs <- cumsum(d$length); k <- which(cs >= q)[1]; if (is.na(k)) NA else round(200 * d$radius[k], 1) }
  qs <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3)
  cat("\nDiametro (cm) alle lunghezze di catena", qs, "\n  porting :", sapply(qs, function(q) qf(r, q)), "\n  originale:", sapply(qs, function(q) qf(o, q)), "\n")
  invisible(list(R = r, O = o))
}
