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
# DBH DIRETTO: misura del diametro a 1.3 m dal suolo senza passare dalla ricostruzione QSM.
# Replica nello spirito il protocollo di misura di riferimento (Pasi): fetta di 20 cm attorno all'altezza,
# trasformata di Hough per cerchi per trovare centro e raggio, poi fit ai minimi quadrati (geometrico, con pesi robusti)
# sui soli punti entro +-4 cm dal cerchio di Hough.
#
# Ipotesi: i punti P (x, y, z) hanno z = quota dal suolo e x, y centrati vicino all'asse del tronco
# (il centro deve stare entro +-10 cm dall'origine: vedi `off` in hough_cerchio; con ricentra = TRUE si centra sulla mediana).
# Parametri fissati prima di guardare i risultati (non tarati sul riferimento di campo).
# NOTA: scritto in Python e verificato li' (21 alberi, mediana errore assoluto 1.3%); questa versione R NON e' stata
# eseguita da chi l'ha scritta: va confrontata con la colonna `atteso_python` di valida_dbh_diretto.R.
# ============================================================

kasa_cerchio <- function(x, y) {
  A <- cbind(x, y, 1); b <- -(x^2 + y^2)
  s <- tryCatch(qr.solve(A, b), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  cx <- -s[1] / 2; cy <- -s[2] / 2; r2 <- cx^2 + cy^2 - s[3]
  if (!is.finite(r2) || r2 <= 0) return(NULL)
  c(cx = unname(cx), cy = unname(cy), r = unname(sqrt(r2)))
}

hough_cerchio <- function(xy, r_min = 0.05, r_max = 0.40, db = 0.005, off = 0.10, passo = 0.01) {
  brk <- seq(r_min, r_max + db, by = db); nb <- length(brk) - 1
  migliore <- -1; out <- NULL
  for (dx in seq(-off, off, by = passo)) {
    for (dy in seq(-off, off, by = passo)) {
      d <- sqrt((xy[, 1] - dx)^2 + (xy[, 2] - dy)^2)
      h <- tabulate(findInterval(d, brk), nbins = nb)
      hs <- h[1:(nb - 2)] + h[2:(nb - 1)] + h[3:nb]            # finestra mobile di 3 classi (1.5 cm)
      k <- which.max(hs)
      if (hs[k] > migliore) { migliore <- hs[k]; out <- c(dx, dy, brk[k + 1] + db / 2) }
    }
  }
  out                                                             # centro x, centro y, raggio
}

cerchio_geometrico <- function(x, y, p0, n_iter = 30, scala = 0.01) {
  p <- p0
  for (it in seq_len(n_iter)) {
    dx <- x - p[1]; dy <- y - p[2]; d <- sqrt(dx^2 + dy^2); d[d < 1e-9] <- 1e-9
    res <- d - p[3]
    J <- cbind(-dx / d, -dy / d, -1)
    w <- 1 / sqrt(1 + (res / scala)^2)                            # pesi robusti (tipo soft-L1)
    JW <- J * w
    step <- tryCatch(qr.solve(crossprod(JW, J), crossprod(JW, res)), error = function(e) NULL)
    if (is.null(step)) break
    p <- p - as.numeric(step)
    if (max(abs(step)) < 1e-6) break
  }
  p
}

dbh_diretto <- function(P, h = 1.3, spessore = 0.2, raggio_max = 0.6, tol = 0.04, ricentra = FALSE) {
  m <- P[, 3] >= h - spessore / 2 & P[, 3] < h + spessore / 2 & sqrt(P[, 1]^2 + P[, 2]^2) < raggio_max
  xy <- P[m, 1:2, drop = FALSE]
  if (nrow(xy) < 50) return(list(diametro_cm = NA_real_, hough_cm = NA_real_, kasa_cm = NA_real_, n = nrow(xy), n_fit = 0))
  c0 <- c(0, 0)
  if (ricentra) { c0 <- c(median(xy[, 1]), median(xy[, 2])); xy <- sweep(xy, 2, c0) }
  ks <- kasa_cerchio(xy[, 1], xy[, 2])
  hz <- hough_cerchio(xy)
  sel <- abs(sqrt((xy[, 1] - hz[1])^2 + (xy[, 2] - hz[2])^2) - hz[3]) < tol
  p <- if (sum(sel) >= 20) cerchio_geometrico(xy[sel, 1], xy[sel, 2], hz) else hz
  list(diametro_cm = 200 * p[3], hough_cm = 200 * hz[3], kasa_cm = if (is.null(ks)) NA_real_ else 200 * ks[["r"]],
       centro = p[1:2] + c0, n = nrow(xy), n_fit = sum(sel))
}
