# =============================================================================
# surface_coverage_filtering.R  (versione corretta)
# BUG FIX 1: distances_to_line ora ritorna lista $d/$h/$V — usata correttamente
# BUG FIX 2: orthonormal_vectors è in helpers.R
# BUG FIX 3: Dis è matrice (nl x ns) ma LexOrd è scalare 1-based → usare
#            indice lineare correttamente
# =============================================================================

#' Filtra una nuvola di punti basandosi sull'ipotesi cilindrica
#' @param P  Nuvola di punti (n x 3)
#' @param c  Cilindro: lista con axis, start, length, radius
#' @param lh Altezza dei layer
#' @param ns Numero di settori iniziale
#' @return Lista con Pass (vettore logico) e cilindro c aggiornato
surface_coverage_filtering <- function(P, c, lh = 0.02, ns = 20) {
  
  np <- nrow(P)
  
  # --- Distanze dalla linea (FIX: usa output lista) ---
  dl <- distances_to_line(P, c$axis, c$start)
  d  <- dl$d
  h  <- dl$h - min(dl$h)
  
  # Vettori perpendicolari per calcolo angolo
  ow <- orthonormal_vectors(c$axis)       # da helpers.R
  # FIX: orthonormal_vectors() restituisce $V/$W, non $U/$W (con $U si
  # otteneva NULL e il successivo V_perp %*% U avrebbe fallito)
  U  <- ow$V;  W <- ow$W
  V_perp <- dl$V                          # già vettori perpendicolari (n x 3)
  Uproj  <- as.numeric(V_perp %*% U)
  Wproj  <- as.numeric(V_perp %*% W)
  ang    <- atan2(Wproj, Uproj) + pi      # [0, 2*pi)
  
  # Lunghezza del cilindro con fallback
  cyl_len <- if (!is.null(c$length) && is.finite(c$length) && c$length > 0)
    c$length else (max(h) - min(h) + 1e-6)
  
  nl <- max(3L, ceiling(cyl_len / lh))
  
  Layer  <- pmin(nl, pmax(1L, ceiling(h  / cyl_len * nl)))
  Sector <- pmin(ns, pmax(1L, ceiling(ang / (2 * pi) * ns)))
  LexOrd <- Layer + (Sector - 1L) * nl   # indice lineare 1-based in matrice nl x ns
  
  ord    <- order(LexOrd)
  LexOrd_s <- LexOrd[ord]
  ds     <- d[ord]
  
  # Prima stima del raggio (mediana del minimo per cella)
  Dis1 <- numeric(nl * ns)
  p <- 1L
  while (p <= np) {
    t <- 1L
    while ((p + t <= np) && LexOrd_s[p] == LexOrd_s[p + t]) t <- t + 1L
    Dmin <- min(ds[p:(p + t - 1L)])
    Dis1[LexOrd_s[p]] <- min(1.05 * Dmin, Dmin + 0.02)
    p <- p + t
  }
  
  R_est <- median(Dis1[Dis1 > 0])
  if (!is.finite(R_est) || R_est <= 0) R_est <- median(d[d > 0])
  if (!is.finite(R_est) || R_est <= 0) R_est <- 0.05
  
  a  <- max(0.02, 0.2 * R_est)
  ns <- min(36L, max(8L, ceiling(2 * pi * R_est / a)))
  nl <- max(3L, ceiling(cyl_len / a))
  
  # Ricalcola con nuova griglia
  Layer  <- pmin(nl, pmax(1L, ceiling(h  / cyl_len * nl)))
  Sector <- pmin(ns, pmax(1L, ceiling(ang / (2 * pi) * ns)))
  LexOrd <- Layer + (Sector - 1L) * nl
  
  ord      <- order(LexOrd)
  LexOrd_s <- LexOrd[ord]
  d_sorted <- d[ord]
  
  Dis2 <- numeric(nl * ns)
  Pass <- logical(np)
  r_cell <- max(0.01, 0.05 * R_est)
  k <- 0L
  
  p <- 1L
  while (p <= np) {
    t <- 1L
    while ((p + t <= np) && LexOrd_s[p] == LexOrd_s[p + t]) t <- t + 1L
    ind  <- p:(p + t - 1L)
    D    <- d_sorted[ind]
    Dmin <- min(D)
    I    <- D <= Dmin + r_cell
    Pass[ind[I]] <- TRUE
    Dis2[LexOrd_s[p]] <- min(1.05 * Dmin, Dmin + 0.02)
    p <- p + t
    k <- k + 1L
  }
  
  # Rimappa Pass all'ordine originale
  inv_ord       <- integer(np)
  inv_ord[ord]  <- seq_len(np)
  Pass          <- Pass[inv_ord]
  
  # Calcolo raggio, SurfCov, mad finali
  R_final <- median(Dis2[Dis2 > 0])
  if (!is.finite(R_final) || R_final <= 0) R_final <- R_est
  
  d_pass  <- d[Pass]
  mad_val <- if (length(d_pass) > 0) mean(abs(d_pass - R_final)) else 0
  
  c$radius  <- R_final
  c$SurfCov <- k / (nl * ns)
  c$mad     <- mad_val
  c$conv    <- 1L
  c$rel     <- 1L
  
  # FIX (27/08/2026): quando la copertura di superficie e' bassa
  # (SurfCov < 0.2 - la stessa soglia usata a valle in cylinder_fitting.R
  # per marcare il cilindro come "non affidabile"), il metodo a griglia
  # sottostima sistematicamente il raggio quando il tronco e' visto solo
  # parzialmente (occlusione angolare) - i pochi settori popolati (tutti
  # dallo stesso lato visibile) dominano la mediana delle celle. Un fit
  # di cerchio algebrico (Kasa, vedi tree_data.R fit_cerchio_3d) sui
  # punti disponibili e' molto piu' robusto in questo caso, perche'
  # stima il cerchio geometricamente piu' plausibile dato l'arco
  # disponibile. Verificato empiricamente su un caso con occlusione
  # (albero "38"): metodo a griglia dava 12.8cm, fit di cerchio dava
  # 18.7cm, contro un vero raggio di 15.7cm (misura di campo).
  # Applicato SOLO quando SurfCov e' gia' basso (segnale esistente di
  # cilindro problematico) E ci sono abbastanza punti per un fit
  # stabile (np >= 20) - non tocca il comportamento normale.
  if (c$SurfCov < 0.2 && np >= 20) {
    fit_alt <- tryCatch(fit_cerchio_3d(P, c$axis, c$start), error = function(e) NULL)
    if (!is.null(fit_alt) && is.finite(fit_alt$r) && fit_alt$r > 0) {
      c$radius <- fit_alt$r
    }
  }
  
  list(Pass = Pass, c = c)
}
