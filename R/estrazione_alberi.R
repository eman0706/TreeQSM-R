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
# ESTRAZIONE DEI SINGOLI ALBERI DAL LAS DI PARCELLA (versione R)
# Passi: (1) quota dal suolo con un DTM grossolano; (2) rilevamento dei tronchi nella fascia 1.1-1.5 m;
#        (3) tracciamento dell'asse inclinato del fusto; (4) selezione dei punti entro R dall'asse; (5) un punto per cubo di 5 mm.
# Scritto a partire dalla versione Python verificata (21 alberi, DBH diretto 1.4% di errore mediano sulla nuvola 15);
# questa versione R NON e' stata eseguita da chi l'ha scritta: va confrontata con verifica_estrazione.R.
# Richiede i pacchetti lidR e dbscan.
# ============================================================
library(lidR); library(dbscan)

carica_las_xyz <- function(percorso) {
  las <- lidR::readLAS(percorso, select = "xyz")
  data.frame(X = las@data$X, Y = las@data$Y, Z = las@data$Z)
}

# Quota relativa al suolo: 2 percentile di Z per cella (default 2 m), celle vuote riempite con la media delle vicine
dtm_celle <- function(X, Y, Z, cs = 2, q = 0.02) {
  ix <- floor((X - min(X)) / cs); iy <- floor((Y - min(Y)) / cs)
  nx <- max(ix) + 1; ny <- max(iy) + 1
  chiave <- ix * ny + iy
  qz <- tapply(Z, chiave, quantile, probs = q, names = FALSE)
  k <- as.integer(names(qz))
  M <- matrix(NA_real_, nx, ny)
  M[cbind(k %/% ny + 1, k %% ny + 1)] <- as.numeric(qz)
  for (it in 1:8) {
    if (!anyNA(M)) break
    Mn <- M
    for (a in seq_len(nx)) for (b in seq_len(ny)) if (is.na(M[a, b])) {
      v <- M[max(1, a - 1):min(nx, a + 1), max(1, b - 1):min(ny, b + 1)]
      if (any(!is.na(v))) Mn[a, b] <- mean(v, na.rm = TRUE)
    }
    M <- Mn
  }
  Z - M[cbind(ix + 1, iy + 1)]
}

# Tronchi candidati: cluster nella fascia di quota (default 1.1-1.5 m), almeno n_min punti, diametro grezzo plausibile
rileva_tronchi <- function(X, Y, zrel, banda = c(1.1, 1.5), eps = 0.15, minPts = 15, n_min = 100, diam_cm = c(8, 70)) {
  s <- zrel >= banda[1] & zrel < banda[2]
  P2 <- cbind(X[s], Y[s])
  db <- dbscan::dbscan(P2, eps = eps, minPts = minPts)
  out <- NULL
  for (k in setdiff(unique(db$cluster), 0)) {
    p <- P2[db$cluster == k, , drop = FALSE]
    if (nrow(p) < n_min) next
    c0 <- c(median(p[, 1]), median(p[, 2]))
    dm <- 2 * median(sqrt((p[, 1] - c0[1])^2 + (p[, 2] - c0[2])^2)) * 100
    if (dm >= diam_cm[1] && dm <= diam_cm[2]) out <- rbind(out, c(x = c0[1], y = c0[2], n = nrow(p), diam_grezzo_cm = dm))
  }
  as.data.frame(out)
}

# Asse del fusto tracciato per fette di 0.25 m, verso l'alto (direzione = 1) o il basso (-1) da quota z_da
traccia_asse <- function(xy, zr, c0, z_da, direzione, passo = 0.25, rtrack = 0.35, maxstep = 0.2, nmin = 30, zlim = 14) {
  z_c <- numeric(0); c_xy <- matrix(numeric(0), 0, 2); cc <- c0; z0 <- z_da
  while (z0 >= 0 && z0 <= zlim) {
    lo <- if (direzione > 0) z0 else z0 - passo; hi <- lo + passo
    m <- zr >= lo & zr < hi & sqrt((xy[, 1] - cc[1])^2 + (xy[, 2] - cc[2])^2) < rtrack
    if (sum(m) < nmin) break
    cn <- c(median(xy[m, 1]), median(xy[m, 2])); dv <- cn - cc; nd <- sqrt(sum(dv^2))
    if (nd > maxstep) cn <- cc + dv / nd * maxstep
    cc <- cn; z_c <- c(z_c, (lo + hi) / 2); c_xy <- rbind(c_xy, cc)
    z0 <- z0 + direzione * passo
  }
  list(z = z_c, xy = c_xy)
}

# Punti dell'albero: entro R dall'asse tracciato a ogni quota; coordinate relative al centro c0 a 1.3 m; z = quota dal suolo
estrai_albero <- function(X, Y, zrel, c0, R, margine = 3, zmax = 14) {
  dx <- X - c0[1]; dy <- Y - c0[2]
  idx <- which(abs(dx) < R + margine & abs(dy) < R + margine & zrel > -0.05 & zrel < zmax)
  xy <- cbind(dx[idx], dy[idx]); zr <- zrel[idx]
  su <- traccia_asse(xy, zr, c(0, 0), 1.3, +1); giu <- traccia_asse(xy, zr, c(0, 0), 1.3, -1)
  zk <- c(giu$z, su$z); ck <- rbind(giu$xy, su$xy)
  if (length(zk) < 3) { zk <- c(0, 14); ck <- matrix(0, 2, 2) }
  o <- order(zk); zk <- zk[o]; ck <- ck[o, , drop = FALSE]
  cx <- approx(zk, ck[, 1], xout = zr, rule = 2)$y; cy <- approx(zk, ck[, 2], xout = zr, rule = 2)$y
  tieni <- sqrt((xy[, 1] - cx)^2 + (xy[, 2] - cy)^2) < R
  P <- cbind(xy[tieni, , drop = FALSE], zr[tieni])
  P[!duplicated(floor(P / 0.005)), , drop = FALSE]                  # un punto per cubo di 5 mm
}

# Raggio di estrazione: min(1.2, 0.45 x distanza dal tronco rilevato piu' vicino)
raggio_estrazione <- function(c0, tronchi) {
  d <- sort(sqrt((tronchi$x - c0[1])^2 + (tronchi$y - c0[2])^2))
  min(1.2, 0.45 * d[2])
}
