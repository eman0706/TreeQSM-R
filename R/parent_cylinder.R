# This file is part of TREEQSM.
# 
# TREEQSM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# TREEQSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

#' Parent Cylinder
#'
#' Porting fedele di parent_cylinder.m (parte di cylinders.m v3.0.0,
#' Pasi Raumonen). Trova il cilindro genitore tra i possibili candidati
#' del segmento genitore, controllando se l'asse del primo cilindro del
#' segmento corrente, esteso, attraversa la superficie di uno dei
#' cilindri vicini nel segmento genitore (risolvendo l'equazione
#' quadratica cilindro-retta). Aggiusta il cilindro cosi' che parta
#' dalla superficie del genitore.
#'
#' NOTA FEDELTA': un ramo specifico dei sei casi simmetrici (quello con
#' x1>0, x2>0, x2>x1) nell'originale usa h(j,1) invece di h(j,2) nella
#' correzione del segno - sembra un refuso nel sorgente MATLAB originale
#' (ovunque altrove il segno di h coincide con quello di x selezionato),
#' ma e' stato replicato esattamente cosi' com'e'.
#'
#' @param SPar segment$ParentSegment (matrice ns x 2)
#' @param SChi segment$ChildSegment (lista)
#' @param CiS Lista: cilindri appartenenti a ciascun segmento (aggiornata
#'   dal chiamante mano a mano che i cilindri vengono prodotti)
#' @param cylinder Struttura globale dei cilindri gia' prodotti (liste di
#'   vettori/matrici: start, axis, radius, length)
#' @param cyl Cilindri appena fittati per il segmento corrente (lista con
#'   radius, length, start, axis, mad, SurfCov, conv, rel - tutti vettori/
#'   matrici allineati per indice)
#' @param si Indice del segmento corrente
#' @return Lista con PC (indice cilindro genitore, o integer(0) se
#'   nessuno), cyl (aggiornato), added (TRUE se e' stato aggiunto un
#'   cilindro extra di raccordo)
#' @export
parent_cylinder <- function(SPar, SChi, CiS, cylinder, cyl, si) {
  
  rad <- cyl$radius
  len <- cyl$length
  sta <- cyl$start
  axe <- cyl$axis
  
  nc <- length(rad)
  added <- FALSE
  
  par_seg <- SPar[si, 1]
  
  if (par_seg > 0) {
    s <- par_seg
    PC <- CiS[[s]]
    if (is.null(PC)) PC <- integer(0)
    
    if (length(PC) > 1) {
      D <- mat_vec_subtraction(-cylinder$start[PC, , drop = FALSE], -sta[1, ])
      d <- rowSums(D * D)
      ord <- order(d)
      if (length(PC) > 3) ord <- ord[1:4]
      pc <- PC[ord]
      ParentFound <- FALSE
    } else if (length(PC) == 1) {
      ParentFound <- TRUE
    } else {
      PC <- integer(0)
      ParentFound <- TRUE
    }
    
    if (!ParentFound) {
      pc0 <- pc
      n <- length(pc)
      x <- matrix(0, n, 2)
      h <- matrix(0, n, 2)
      Axe <- cylinder$axis[pc, , drop = FALSE]
      Sta <- cylinder$start[pc, , drop = FALSE]
      
      for (j in seq_len(n)) {
        axe1  <- axe[1, ]; sta1 <- sta[1, ]
        Axe_j <- Axe[j, ]; Sta_j <- Sta[j, ]
        A_v <- axe1 - as.numeric(axe1 %*% Axe_j) * Axe_j
        B_v <- sta1 - Sta_j - as.numeric(sta1 %*% Axe_j) * Axe_j +
          as.numeric(Sta_j %*% Axe_j) * Axe_j
        e_  <- sum(A_v * A_v)
        f_  <- 2 * sum(A_v * B_v)
        g_  <- sum(B_v * B_v) - cylinder$radius[pc[j]]^2
        disc <- f_^2 - 4 * e_ * g_
        if (is.finite(disc) && disc >= 0 && abs(e_) > 1e-15) {
          di <- sqrt(disc)
          s1 <- (-f_ + di) / (2 * e_)
          s2 <- (-f_ - di) / (2 * e_)
          x[j, ] <- c(s1, s2)
          h[j, 1] <- sum(sta1 * Axe_j) + s1 * sum(axe1 * Axe_j) - sum(Sta_j * Axe_j)
          h[j, 2] <- sum(sta1 * Axe_j) + s2 * sum(axe1 * Axe_j) - sum(Sta_j * Axe_j)
        }
      }
      
      I <- x[, 1] != 0
      pc <- pc0[I]; x <- x[I, , drop = FALSE]; h <- h[I, , drop = FALSE]
      n <- sum(I)
      X <- matrix(0, n, 3)
      Len <- cylinder$length[pc]
      
      j <- 1L
      while (j <= n && !ParentFound) {
        x1 <- x[j, 1]; x2 <- x[j, 2]; h1 <- h[j, 1]; h2 <- h[j, 2]; Lj <- Len[j]
        
        matched <- TRUE
        if (x1 > 0 && x2 < 0) {
          xs <- x1; hs <- h1
        } else if (x1 < 0 && x2 > 0 && len[1] - x2 > 0) {
          xs <- x2; hs <- h2
        } else if (x1 < 0 && x2 < 0 && x2 < x1 && len[1] - x1 > 0) {
          xs <- x1; hs <- h1
        } else if (x1 < 0 && x2 < 0 && x2 > x1 && len[1] - x2 > 0) {
          xs <- x2; hs <- h2
        } else if (x1 > 0 && x2 > 0 && x2 < x1 && len[1] - x1 > 0) {
          xs <- x1; hs <- h1
        } else if (x1 > 0 && x2 > 0 && x2 > x1 && len[1] - x2 > 0) {
          xs <- x2; hs <- h2
        } else {
          matched <- FALSE
        }
        
        if (matched) {
          if (hs >= 0 && hs <= Lj && len[1] - xs > 0) {
            PC <- pc[j]
            sta[1, ] <- sta[1, ] + xs * axe[1, ]
            len[1] <- len[1] - xs
            ParentFound <- TRUE
          } else if (len[1] - xs > 0) {
            # NOTA FEDELTA': nel ramo x1>0,x2>0,x2>x1 l'originale usa h1
            # invece di h2 qui (vedi nota sopra) - replicato cosi' com'e'
            h_sign <- if (x1 > 0 && x2 > 0 && x2 > x1) h1 else hs
            if (h_sign < 0) X[j, ] <- c(xs, abs(hs), 0) else X[j, ] <- c(xs, hs - Lj, 0)
          } else {
            X[j, ] <- c(xs, hs, 1)
          }
        }
        j <- j + 1L
      }
      
      if (!ParentFound && n > 0) {
        best <- which.min(X[, 2])
        H <- X[best, 2]
        Xrow <- X[best, ]
        if (Xrow[3] == 0 && H < 0.1 * Len[best]) {
          PC <- pc[best]
          sta[1, ] <- sta[1, ] + Xrow[1] * axe[1, ]
          len[1] <- len[1] - Xrow[1]
          ParentFound <- TRUE
        } else {
          PC <- pc[best]
          if (nc > 1 && Xrow[1] <= rad[1] && abs(Xrow[2]) <= 1.25 * cylinder$length[PC]) {
            S <- sta[1, ] + Xrow[1] * axe[1, ]
            V <- sta[2, ] + len[2] * axe[2, ] - S
            Vn <- sqrt(sum(V * V))
            len[2] <- Vn; len <- len[2:nc]
            axe[2, ] <- V / Vn; axe <- axe[2:nc, , drop = FALSE]
            sta[2, ] <- S; sta <- sta[2:nc, , drop = FALSE]
            rad <- rad[2:nc]
            cyl$mad <- cyl$mad[2:nc]
            cyl$SurfCov <- cyl$SurfCov[2:nc]
            nc <- nc - 1L
            ParentFound <- TRUE
          } else if (nc > 1) {
            sta <- sta[2:nc, , drop = FALSE]; len <- len[2:nc]
            axe <- axe[2:nc, , drop = FALSE]; rad <- rad[2:nc]
            cyl$mad <- cyl$mad[2:nc]
            cyl$SurfCov <- cyl$SurfCov[2:nc]
            nc <- nc - 1L
          } else if (length(SChi[[si]]) == 0) {
            nc <- 0L; PC <- integer(0); ParentFound <- TRUE; rad <- numeric(0)
          } else if (Xrow[1] <= rad[1] && abs(Xrow[2]) <= 1.5 * cylinder$length[PC]) {
            sta[1, ] <- sta[1, ] + Xrow[1] * axe[1, ]
            len[1] <- abs(Xrow[1])
            ParentFound <- TRUE
          }
        }
      }
      
      if (!ParentFound) {
        pc <- pc0
        dbl <- distances_between_lines(sta[1, ], axe[1, ],
                                       cylinder$start[pc, , drop = FALSE],
                                       cylinder$axis[pc, , drop = FALSE])
        Dist <- dbl$Dist
        DistOnLines <- dbl$DistOnLines
        
        I2 <- DistOnLines >= 0
        J2 <- DistOnLines <= cylinder$length[pc]
        I2 <- I2 & J2
        if (!any(I2)) {
          I2 <- DistOnLines >= -0.2 * cylinder$length[pc]
          J2 <- DistOnLines <= 1.2 * cylinder$length[pc]
          I2 <- I2 & J2
        }
        if (any(I2)) {
          pc2 <- pc[I2]; Dist2 <- Dist[I2]; DistOnLines2 <- DistOnLines[I2]
          best <- which.min(Dist2)
          DoL <- DistOnLines2[best]; PC <- pc2[best]
          Q <- cylinder$start[PC, ] + DoL * cylinder$axis[PC, ]
          V <- sta[1, ] - Q; L <- sqrt(sum(V * V)); V <- V / L
          a_ang <- acos(max(-1, min(1, sum(V * cylinder$axis[PC, ]))))
          hgt <- sin(a_ang) * L
          if (hgt > 1e-9) {
            S <- Q + cylinder$radius[PC] / hgt * L * V
            L2 <- (hgt - cylinder$radius[PC]) / hgt * L
            if (is.finite(L2) && L2 > 0.01 && L2 / len[1] > 0.2) {
              nc <- nc + 1L
              sta <- rbind(S, sta); rad <- c(rad[1], rad)
              axe <- rbind(V, axe); len <- c(L2, len)
              cyl$mad <- c(cyl$mad[1], cyl$mad)
              cyl$SurfCov <- c(cyl$SurfCov[1], cyl$SurfCov)
              cyl$rel <- c(cyl$rel[1], cyl$rel)
              cyl$conv <- c(cyl$conv[1], cyl$conv)
              added <- TRUE
            }
          }
        } else {
          V <- -mat_vec_subtraction(cylinder$start[pc, , drop = FALSE], sta[1, ])
          L0 <- sqrt(rowSums(V * V))
          V <- V / L0
          Adot <- as.numeric(V %*% axe[1, ])
          best <- which.max(Adot)
          L1 <- L0[best]; PC <- pc[best]; Vb <- V[best, ]
          a_ang <- acos(max(-1, min(1, sum(Vb * cylinder$axis[PC, ]))))
          hgt <- sin(a_ang) * L1
          if (hgt > 1e-9) {
            S <- cylinder$start[PC, ] + cylinder$radius[PC] / hgt * L1 * Vb
            L2 <- (hgt - cylinder$radius[PC]) / hgt * L1
            if (is.finite(L2) && L2 > 0.01 && L2 / len[1] > 0.2) {
              nc <- nc + 1L
              sta <- rbind(S, sta); rad <- c(rad[1], rad)
              axe <- rbind(Vb, axe); len <- c(L2, len)
              cyl$mad <- c(cyl$mad[1], cyl$mad)
              cyl$SurfCov <- c(cyl$SurfCov[1], cyl$SurfCov)
              cyl$rel <- c(cyl$rel[1], cyl$rel)
              cyl$conv <- c(cyl$conv[1], cyl$conv)
              added <- TRUE
            }
          }
        }
      }
    }
  } else {
    PC <- integer(0)
  }
  
  if (nc > 0) {
    cyl$radius  <- rad[seq_len(nc)]
    cyl$length  <- len[seq_len(nc)]
    cyl$start   <- sta[seq_len(nc), , drop = FALSE]
    cyl$axis    <- axe[seq_len(nc), , drop = FALSE]
    cyl$mad     <- cyl$mad[seq_len(nc)]
    cyl$SurfCov <- cyl$SurfCov[seq_len(nc)]
    cyl$conv    <- cyl$conv[seq_len(nc)]
    cyl$rel     <- cyl$rel[seq_len(nc)]
  } else {
    cyl$radius <- numeric(0); cyl$length <- numeric(0)
    cyl$start <- matrix(0, 0, 3); cyl$axis <- matrix(0, 0, 3)
    cyl$mad <- numeric(0); cyl$SurfCov <- numeric(0)
    cyl$conv <- logical(0); cyl$rel <- logical(0)
  }
  
  list(PC = PC, cyl = cyl, added = added)
}