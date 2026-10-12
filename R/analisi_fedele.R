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

# ANALISI DELLE CORSE FEDELI: DBHcyl del QSM (base fedele + Gauss-Newton), tre seed, entrambe le nuvole, contro il campo, contro la procedura di Greco et al.
# e contro il porting di produzione precedente (modelli_quote_vol.rds). Nessuna casualita': solo lettura dei file; il bootstrap usa set.seed(1).
# PREVISIONI SCRITTE PRIMA: P_R1 mediana |errore| del DBHcyl fedele tra 3.0% e 4.2% su entrambe le nuvole; P_R2 il QSM fedele non e' piu' accurato della procedura
#  di Greco et al. (mediana |err| 3.00% nuvola 15, 3.07% nuvola 11) su nessuna delle due; P_R3 il Wilcoxon appaiato fedele contro produzione precedente non e' significativo (p > 0.05) su entrambe;
#  P_R4 il volume QSM del tratto 1.3-4.0 m resta sottostimato (segno mediano negativo) su entrambe, con mediana |err| entro 3 punti dalla precedente (nuvola 15: 8-12%).
# Uso:  source(".../analisi_fedele.R");  a <- analisi_fedele(c15, c11)      (dopo r15 <- raccolta_modelli_quote(...fedele.rds) e r11, oppure leggendo i file)
analisi_fedele <- function(c15, c11, f15 = file.path(c15, "modelli_quote_fedele.rds"), f11 = file.path(c11, "modelli_quote_fedele.rds"), B = 10000) {
  set.seed(1); out <- list()
  for (nv in c("15", "11")) {
    cart <- if (nv == "15") c15 else c11; fed <- readRDS(if (nv == "15") f15 else f11)$modelli
    vec <- readRDS(file.path(cart, "modelli_quote_vol.rds"))$modelli
    rif <- read.csv(file.path(cart, "riferimento_alberi.csv"), stringsAsFactors = FALSE); rif <- rif[rif$fonte_ref == "campo", ]
    au <- read.csv(file.path(cart, if (nv == "15") "volumi_autori_las15.csv" else "volumi_autori_las11.csv"), stringsAsFactors = FALSE)
    pa <- function(M, n) { d <- M[M$N == n & is.finite(M$dbh_matlab_cm) & M$dbh_matlab_cm >= 5, ]; if (nrow(d) == 0) NA_real_ else median(tapply(d$dbh_matlab_cm, d$off, median)) }
    tb <- data.frame(N = rif$N, ref = rif$ref_cm, fedele = sapply(rif$N, function(n) pa(fed, n)), vecchio = sapply(rif$N, function(n) pa(vec, n)), greco = if ("DiamCyl130_cm" %in% names(au)) au$DiamCyl130_cm[match(rif$N, au$TreeID)] else rep(NA_real_, nrow(rif)))
    for (k in c("fedele", "vecchio", "greco")) tb[[paste0("e_", k)]] <- 100 * (tb[[k]] - tb$ref) / tb$ref
    cat(sprintf("\n=================== NUVOLA %s (DBHcyl, cm; errori %%) ===================\n", nv)); print(round(tb, 1), row.names = FALSE)
    st <- function(e) c(med_abs = median(abs(e), na.rm = TRUE), media_abs = mean(abs(e), na.rm = TRUE), segno = median(e, na.rm = TRUE), oltre15 = sum(abs(e) > 15, na.rm = TRUE), n = sum(is.finite(e)))
    cat("\n"); print(round(rbind(fedele = st(tb$e_fedele), produzione_precedente = st(tb$e_vecchio), Greco = st(tb$e_greco)), 2))
    w <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 6) return(NA_real_); suppressWarnings(wilcox.test(abs(a[ok]), abs(b[ok]), paired = TRUE, exact = FALSE)$p.value) }
    cat(sprintf("Wilcoxon appaiato sugli |errori|: fedele-Greco p = %.3f | fedele-precedente p = %.3f\n", w(tb$e_fedele, tb$e_greco), w(tb$e_fedele, tb$e_vecchio)))
    if (!("DiamCyl130_cm" %in% names(au))) cat(sprintf("  NOTA: il file dei volumi di Greco della nuvola %s non ha la colonna DiamCyl130_cm (colonne: %s): confronto per albero con Greco omesso\n", nv, paste(names(au), collapse = ", ")))
    ok <- is.finite(tb$e_fedele) & is.finite(tb$e_greco)
    if (sum(ok) >= 6) { bs <- replicate(B, { i <- sample(which(ok), replace = TRUE); c(median(abs(tb$e_fedele[i])), median(abs(tb$e_fedele[i])) - median(abs(tb$e_greco[i]))) })
    cat(sprintf("Bootstrap (%d, sugli alberi): mediana |err| fedele %.2f%% [IC95 %.2f; %.2f] | differenza fedele - Greco %+.2f punti [%.2f; %.2f]\n", B, st(tb$e_fedele)[["med_abs"]], quantile(bs[1, ], .025), quantile(bs[1, ], .975), median(bs[2, ]), quantile(bs[2, ], .025), quantile(bs[2, ], .975)))
    }
    a1 <- analizza_quote(list(modelli = fed), cart, NULL, nv); a2 <- analizza_volumi(list(modelli = fed), cart, file.path(cart, if (nv == "15") "volumi_autori_las15.csv" else "volumi_autori_las11.csv"), nv)
    cat("  --- stessi volumi con la corsa PRECEDENTE (per confronto) ---\n"); invisible(analizza_volumi(list(modelli = vec), cart, file.path(cart, if (nv == "15") "volumi_autori_las15.csv" else "volumi_autori_las11.csv"), nv))
    out[[nv]] <- list(tabella = tb, quote = a1, volumi = a2)
  }
  invisible(out)
}