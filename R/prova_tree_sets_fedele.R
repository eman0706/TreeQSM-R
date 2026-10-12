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
# PROVA: porting fedele delle tre funzioni di tree_sets che deviano da MATLAB
#
# diagnosi_tree_sets.R ha mostrato che nel secondo passaggio il grafo dei cover set e' molto frammentato (500-1200 componenti,
# 55-80% dei cover set nella piu' grande) e che define_trunk, define_main_branches e make_tree_connected non aggiungono nessun
# collegamento: nel porting define_main_branches marca come Trunk anche i frammenti staccati e non ripara nulla (in MATLAB ricollega
# le componenti di ogni ramo principale e unisce i rami di primo ordine al fusto); la base del secondo passaggio ha 5 cover set
# invece di 66-147. Qui le tre funzioni sono riscritte da MATLAB (define_base_forb con OnlyTree, define_trunk,
# define_main_branches) e sostituite TEMPORANEAMENTE; si confrontano con la corsa v4 (stessi alberi, stessi seed).
#   configurazione "mb":   solo define_main_branches fedele
#   configurazione "tutte": define_base_forb (OnlyTree), define_trunk e define_main_branches fedeli
# Non si modifica nessun file; tutto viene ripristinato alla fine. NOTA: scritto senza poterlo eseguire; se un passo da' errore
# lo scrive e si ferma.
#
# PREVISIONI SCRITTE PRIMA
#  P_F1: con define_main_branches fedele la quota di cover set nella componente connessa piu' grande (secondo passaggio) sale a >= 95%.
#  P_F2: il DBH (DBHcyl fedele, mediana dei modelli e dei seed) cambia al massimo di 0.7 punti (mediana sui 5 alberi) rispetto
#        alla corsa v4: la riparazione non tocca il fusto vicino alla base.
#  P_F3: l'altezza del QSM (TreeHeight, mediana dei modelli) aumenta di almeno il 5% in almeno 3 alberi su 5.
#  P_F4 (esplorativa): con le tre funzioni fedeli ("tutte") la mediana |errore| di DBHcyl sui 5 alberi resta <= 4.0%.
#
# Uso (sessione ripristinata, surface_coverage corretta):
#   source(".../prova_tree_sets_fedele.R");  r <- prova_tree_sets_fedele(".../alberi_corretti")      # circa 40-60 minuti
# ============================================================
pdist2_r <- function(A, B) { d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B); sqrt(pmax(d2, 0)) }

define_base_forb_fedele <- function(P, cover, aux, inputs, segment = NULL) {
  Ce <- aux$Ce
  rng2 <- function(M) max(apply(M, 2, function(x) diff(range(x))))                 # MATLAB: max(max(.)-min(.)) su x e y
  if (is.null(segment)) {                                                          # prima chiamata (nargin == 4)
    BaseHeight <- min(1.5, 0.02 * aux$Height)
    Base <- aux$Ind[Ce[, 3] < aux$Hmin + BaseHeight]
    if (length(Base) == 0) Base <- aux$Ind[which.min(Ce[, 3])]
    Forb <- aux$Fal
    Wb <- rng2(Ce[Base, 1:2, drop = FALSE]); Wt <- rng2(Ce[, 1:2, drop = FALSE]); k <- 1L
    while (k <= 5L && Wb > 0.3 * Wt) {
      BaseHeight <- max(BaseHeight - 0.05, 0.05)
      Base <- aux$Ind[Ce[, 3] < aux$Hmin + BaseHeight]
      if (length(Base) == 0) Base <- aux$Ind[which.min(Ce[, 3])]
      Wb <- rng2(Ce[Base, 1:2, drop = FALSE]); k <- k + 1L
    }
    return(list(Base = Base, Forb = Forb, cover = cover))
  }
  BaseHeight <- min(1.5, 0.02 * aux$Height)                                        # seconda chiamata (OnlyTree con segment)
  SoP <- segment$SegmentOfPoint[cover$center]
  stem <- aux$Ind[!is.na(SoP) & SoP == 1]
  Base <- stem[Ce[stem, 3] < aux$Hmin + BaseHeight]
  if (length(Base) == 0) Base <- if (length(stem) > 0) stem[which.min(Ce[stem, 3])] else aux$Ind[which.min(Ce[, 3])]
  list(Base = Base, Forb = aux$Fal, cover = cover)
}

define_trunk_fedele <- function(cover, aux, Base, Forb, inputs) {
  Nei <- cover$neighbor; Ce <- aux$Ce; nb <- aux$nb
  Trunk <- aux$Fal; Trunk[Base] <- TRUE
  Exp <- unique(c(Base, unlist(Nei[Base], use.names = FALSE)))
  Exp <- Exp[!Trunk[Exp] | !Forb[Exp]]
  Trunk[Exp] <- TRUE
  L <- 0.25
  H <- max(Ce[Trunk, 3]) - L
  FirstMod <- TRUE; giri <- 0L
  while (length(Exp) > 0 && giri < 20L * nb) {
    giri <- giri + 1L
    H0 <- H; Exp0 <- Exp
    Exp <- unique(c(Exp, unlist(Nei[Exp], use.names = FALSE)))
    Exp <- Exp[!Trunk[Exp]]
    Exp <- Exp[Ce[Exp, 3] >= H]
    Trunk[Exp] <- TRUE
    if (length(Exp) > 0) H <- max(Ce[Exp, 3]) - L
    if ((length(Exp) == 0 || H < H0 + inputs$PatchDiam1 / 2) && H < aux$Height - 5) {
      if (FirstMod) {
        FirstMod <- FALSE
        Min <- apply(Ce[, 1:2, drop = FALSE], 2, min); Max <- apply(Ce[, 1:2, drop = FALSE], 2, max)
        EL <- 0.2; NRect <- ceiling((Max - Min) / EL) + 1
        px <- floor((Ce[, 1] - Min[1]) / EL) + 1; py <- floor((Ce[, 2] - Min[2]) / EL) + 1
      }
      Region <- if (length(Exp) > 0) Exp else Exp0
      X1 <- min(px[Region]); if (X1 <= 2) X1 <- 3
      X2 <- max(px[Region]); if (X2 >= NRect[1] - 1) X2 <- NRect[1] - 2
      Y1 <- min(py[Region]); if (Y1 <= 2) Y1 <- 3
      Y2 <- max(py[Region]); if (Y2 >= NRect[2] - 1) Y2 <- NRect[2] - 2
      K <- px >= X1 - 2 & px <= X2 + 2 & py >= Y1 - 2 & py <= Y2 + 2
      I <- Ce[, 3] > H & Ce[, 3] < H + 2 & K; I[Trunk] <- FALSE
      SetsAbove <- which(I)
      if (length(SetsAbove) > 0) {
        n <- length(Region); m <- length(SetsAbove)
        Dist <- matrix(0, n, m); Cos <- matrix(0, n, m)
        for (i in seq_len(n)) {
          V <- sweep(Ce[SetsAbove, , drop = FALSE], 2, Ce[Region[i], ])
          Len <- rowSums(V * V)
          Dist[i, ] <- Len; Cos[i, ] <- V[, 3] / sqrt(pmax(Len, 1e-300))
        }
        Isel <- Cos > 0.7; t <- 0
        while (!any(Isel) && t < 20) { t <- t + 1; Isel <- Cos > 0.7 - t * 0.05 }
        Dist[!Isel] <- 3
        if (n > 1 && m > 1) {
          Iv <- apply(Dist, 2, which.min); d <- apply(Dist, 2, min); J <- which.min(d); I1 <- Iv[J]
        } else if (n == 1 && m > 1) { J <- which.min(Dist[1, ]); I1 <- 1 } else { I1 <- 1; J <- 1 }
        a <- Region[I1]; b <- SetsAbove[J]
        Nei[[a]] <- c(Nei[[a]], b); Nei[[b]] <- c(Nei[[b]], a)
        Exp <- unique(c(Region, unlist(Nei[Region], use.names = FALSE)))
        Exp <- Exp[!Trunk[Exp]]; Exp <- Exp[Ce[Exp, 3] >= H]
        Trunk[Exp] <- TRUE
        if (length(Exp) > 0) H <- max(Ce[Exp, 3]) - L
      }
    }
  }
  cover$neighbor <- Nei
  list(Trunk = Trunk, cover = cover)
}

define_main_branches_fedele <- function(cover, segment, aux, inputs) {
  Bal <- cover$ball; Nei <- cover$neighbor; Ce <- aux$Ce; nb <- length(Bal)
  SoP <- segment$SegmentOfPoint; maxseg <- max(SoP, na.rm = TRUE)
  if (!is.finite(maxseg) || maxseg < 1) return(list(Trunk = aux$Fal, cover = cover))
  MBI <- logical(maxseg); MBI[1] <- TRUE
  for (idx in c(segment$branch1indexes, segment$branch2indexes, segment$branch3indexes))
    if (!is.na(idx) && idx >= 1 && idx <= maxseg) MBI[idx] <- TRUE
  MainBranches <- integer(nb)
  for (i in seq_len(nb)) {
    bi <- SoP[Bal[[i]]]; bi <- bi[!is.na(bi) & bi > 0]
    if (length(bi) > 0) { ind <- min(bi); if (MBI[ind]) MainBranches[i] <- ind }
  }
  Trunk <- aux$Fal; Trunk[MainBranches > 0] <- TRUE
  cp <- cubical_partition(Ce, 3 * inputs$PatchDiam2Max, 10); Par <- cp$Partition; CC <- cp$CubeCoord; N <- dim(Par)
  BI <- max(MainBranches)
  for (i in seq_len(BI)) if (MBI[i]) {
    Branch <- MainBranches == i
    Comps <- connected_components(Nei, Branch, 1, aux$Fal)$Components; n <- length(Comps); ciclo <- 0L
    while (n > 1 && ciclo < 100L) {
      ciclo <- ciclo + 1L
      for (j in seq_len(n)) {
        comp <- Comps[[j]]
        cu <- unique(CC[comp, , drop = FALSE]); m <- nrow(cu)
        t <- 0L; NearSets <- integer(0)
        while (length(NearSets) == 0 && t < max(N)) {
          t <- t + 1L; flag <- aux$Fal
          for (k in seq_len(m)) {
            x1 <- max(1, cu[k, 1] - t); x2 <- min(cu[k, 1] + t, N[1]); y1 <- max(1, cu[k, 2] - t); y2 <- min(cu[k, 2] + t, N[2])
            z1 <- max(1, cu[k, 3] - t); z2 <- min(cu[k, 3] + t, N[3])
            balls <- unlist(Par[x1:x2, y1:y2, z1:z2], use.names = FALSE)
            balls <- balls[Branch[balls]]
            flag[balls] <- TRUE
          }
          flag[comp] <- FALSE
          NearSets <- which(flag)
        }
        if (length(NearSets) > 0) {
          d <- pdist2_r(Ce[comp, , drop = FALSE], Ce[NearSets, , drop = FALSE])
          ij <- which(d == min(d), arr.ind = TRUE)[1, ]
          a <- comp[ij[1]]; b <- NearSets[ij[2]]
          Nei[[a]] <- c(Nei[[a]], b); Nei[[b]] <- c(Nei[[b]], a)
        }
      }
      Comps <- connected_components(Nei, Branch, 1, aux$Fal)$Components; n <- length(Comps)
    }
  }
  Stem <- which(MainBranches == 1)
  b1 <- segment$branch1indexes; b1 <- b1[!is.na(b1) & b1 >= 1 & b1 <= maxseg]
  BI2 <- if (length(b1) > 0) max(b1) else 0
  MBI2 <- logical(maxseg); MBI2[b1] <- TRUE
  if (BI2 >= 2 && length(Stem) > 0) for (i in 2:BI2) if (MBI2[i]) {
    Branch <- which(MainBranches == i)
    if (length(Branch) > 0) {
      vic <- unlist(Nei[Branch], use.names = FALSE)
      if (!any(MainBranches[vic] == 1)) {
        d <- pdist2_r(Ce[Branch, , drop = FALSE], Ce[Stem, , drop = FALSE])
        ij <- which(d == min(d), arr.ind = TRUE)[1, ]
        a <- Branch[ij[1]]; b <- Stem[ij[2]]
        Nei[[a]] <- c(Nei[[a]], b); Nei[[b]] <- c(Nei[[b]], a)
      }
    }
  }
  cover$neighbor <- Nei
  cc <- connected_components(cover$neighbor, Trunk, 0, aux$Fal); comps <- cc$Components; cs <- cc$CompSize
  if (length(cs) > 1) {
    o <- order(cs, decreasing = TRUE); comps <- comps[o]; cs <- cs[o]
    StemV <- MainBranches == 1; Trunk <- aux$Fal; i <- 1L; C <- comps[[i]]
    while (i < length(cs) && !any(StemV[C])) { i <- i + 1L; C <- comps[[i]] }
    Trunk[C] <- TRUE
  }
  list(Trunk = Trunk, cover = cover)
}

prova_tree_sets_fedele <- function(cartella, cartella_R = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", alberi = c(4, 35, 1, 16, 28), offsets = 0:2, config = c("mb", "tutte"),
                                   file_baseline = file.path(cartella, "modelli_v4_sc_corretta.rds")) {
  stopifnot(exists("treeqsm"), exists("dbh_robusto"), exists("inputs_prod"), exists("dbh_cylinder"), exists("define_trunk"), exists("define_main_branches"),
            exists("define_base_forb"), exists("cubical_partition"), exists("connected_components"))
  if (!("full" %in% names(formals(surface_coverage)))) stop("surface_coverage non e' la versione corretta")
  orig <- list(define_trunk = get("define_trunk", envir = globalenv()), define_main_branches = get("define_main_branches", envir = globalenv()),
               define_base_forb = get("define_base_forb", envir = globalenv()), treeqsm = get("treeqsm", envir = globalenv()),
               relative_size = get("relative_size", envir = globalenv()))
  on.exit(for (n in names(orig)) assign(n, orig[[n]], envir = globalenv()), add = TRUE)
  source(file.path(cartella_R, "relative_size.R"))
  txt <- readLines(file.path(cartella_R, "treeqsm.R"), warn = FALSE)
  i_td <- grep("treedata = tdata\\$treedata,", txt); if (length(i_td) != 1) stop("treeqsm.R: punto di inserimento non trovato")
  txt[i_td] <- paste0(txt[i_td], "\n            trunk_points = trunk_points,")
  tf <- tempfile(fileext = ".R"); writeLines(txt, tf); source(tf)
  assign(".off", 0L, envir = globalenv())
  assign("set.seed", function(seed, ...) base::set.seed(as.integer((as.numeric(seed) + get(".off", envir = globalenv()) * 7919) %% 2147483647), ...), envir = globalenv())
  on.exit(if (exists("set.seed", envir = globalenv(), inherits = FALSE)) rm("set.seed", envir = globalenv()), add = TRUE)
  
  comp <- function(Nei) { n <- length(Nei); lab <- integer(n); k <- 0L
  for (s in seq_len(n)) if (lab[s] == 0L) { k <- k + 1L; lab[s] <- k; pila <- s
  while (length(pila)) { u <- pila[length(pila)]; pila <- pila[-length(pila)]; v <- Nei[[u]]; v <- v[!is.na(v) & v >= 1 & v <= n]; v <- unique(v[lab[v] == 0L])
  if (length(v)) { lab[v] <- k; pila <- c(pila, v) } } }
  tab <- tabulate(lab, k); c(k = k, pct = 100 * max(tab) / n) }
  frazioni <- new.env(); frazioni$v <- numeric(0)
  gv <- function(x) if (is.null(x) || length(x) == 0 || !is.finite(x[1])) NA_real_ else as.numeric(x[1])
  dbh_fed <- function(m) { if (is.null(m$trunk_points)) return(NA_real_); cy <- m$cylinder
  r <- tryCatch(dbh_cylinder(cy$branch == 1, cy, seq_along(cy$radius), m$trunk_points), error = function(e) NULL); if (is.null(r)) NA_real_ else r$DBHcyl }
  rif <- read.csv(file.path(cartella, "riferimento_alberi.csv"), stringsAsFactors = FALSE); v <- setNames(rif$ref_cm, rif$N)
  base <- readRDS(file_baseline)
  
  configura <- function(cfg) {
    assign("define_trunk", orig$define_trunk, envir = globalenv()); assign("define_base_forb", orig$define_base_forb, envir = globalenv())
    mb <- if (cfg %in% c("mb", "mb+base", "mb+trunk", "tutte")) define_main_branches_fedele else orig$define_main_branches
    assign("define_main_branches", function(cover, segment, aux, inputs) {
      out <- mb(cover, segment, aux, inputs); frazioni$v <- c(frazioni$v, comp(out$cover$neighbor)[["pct"]]); out }, envir = globalenv())
    if (cfg %in% c("tutte", "mb+trunk")) assign("define_trunk", define_trunk_fedele, envir = globalenv())
    if (cfg %in% c("tutte", "mb+base")) assign("define_base_forb", define_base_forb_fedele, envir = globalenv())
  }
  risultati <- list()
  for (cfg in config) {
    configura(cfg)
    for (n in alberi) {
      P <- as.matrix(read.csv(file.path(cartella, sprintf("albero_%03d.csv", n)))); mode(P) <- "numeric"; Pp <- P[P[, 3] >= 0.15, , drop = FALSE]
      frazioni$v <- numeric(0); t0 <- Sys.time(); Fs <- numeric(0); As <- numeric(0); Hs <- numeric(0); Cs <- numeric(0)
      for (o in offsets) {
        assign(".off", as.integer(o), envir = globalenv())
        q <- NULL; msg <- NA_character_
        invisible(capture.output(q <- tryCatch(treeqsm(Pp, inputs_prod), error = function(e) { msg <<- conditionMessage(e); NULL })))
        if (is.null(q) || length(q) == 0) { cat("  errore treeqsm:", msg, "\n"); next }
        As <- c(As, tryCatch(dbh_robusto(q), error = function(e) NA_real_))
        f <- 100 * vapply(q, function(m) gv(dbh_fed(m)), numeric(1)); Fs <- c(Fs, if (any(f >= 5, na.rm = TRUE)) median(f[f >= 5], na.rm = TRUE) else NA_real_)
        Hs <- c(Hs, median(vapply(q, function(m) gv(if (!is.null(m$treedata$TreeHeight)) m$treedata$TreeHeight else m$treedata$Height), numeric(1)), na.rm = TRUE))
        Cs <- c(Cs, median(vapply(q, function(m) length(m$cylinder$radius), numeric(1))))
      }
      mb0 <- base$modelli[base$modelli$N == n, ]
      risultati[[length(risultati) + 1]] <- data.frame(cfg = cfg, N = n, ref_cm = unname(v[as.character(n)]),
                                                       F_v4 = round(median(tapply(mb0$dbh_matlab_cm[mb0$dbh_matlab_cm >= 5], mb0$off[mb0$dbh_matlab_cm >= 5], median)), 1), F_nuovo = round(median(Fs, na.rm = TRUE), 1),
                                                       h_v4 = round(median(mb0$altezza, na.rm = TRUE), 2), h_nuova = round(median(Hs, na.rm = TRUE), 2),
                                                       cil_v4 = median(mb0$n_cil), cil_nuovi = median(Cs), comp_grande_pct = round(median(frazioni$v), 1), minuti = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
      cat(sprintf("[%s] %s albero %d fatto\n", format(Sys.time(), "%H:%M:%S"), cfg, n))
    }
  }
  tb <- do.call(rbind, risultati)
  tb$errF_v4 <- round(100 * (tb$F_v4 - tb$ref_cm) / tb$ref_cm, 1); tb$errF_nuovo <- round(100 * (tb$F_nuovo - tb$ref_cm) / tb$ref_cm, 1)
  tb$var_h_pct <- round(100 * (tb$h_nuova - tb$h_v4) / tb$h_v4, 1)
  print(tb[, c("cfg", "N", "ref_cm", "F_v4", "F_nuovo", "errF_v4", "errF_nuovo", "h_v4", "h_nuova", "var_h_pct", "cil_v4", "cil_nuovi", "comp_grande_pct", "minuti")], row.names = FALSE)
  cat("\nConfigurazioni disponibili: mb (solo define_main_branches), mb+base (+ define_base_forb), mb+trunk (+ define_trunk), tutte (le tre)\n")
  for (cf in unique(tb$cfg)) { g <- tb[tb$cfg == cf, ]
  cat(sprintf("%-9s componente grande mediana %.1f%% | variazione mediana DBHcyl %.1f punti (massima %.1f) | altezza QSM +%.1f%% (mediana) | mediana |errore| DBHcyl %.1f%% (v4: %.1f%%)\n",
              cf, median(g$comp_grande_pct), median(abs(g$errF_nuovo - g$errF_v4)), max(abs(g$errF_nuovo - g$errF_v4)), median(g$var_h_pct), median(abs(g$errF_nuovo)), median(abs(g$errF_v4)))) }
  if (all(c("mb", "tutte") %in% tb$cfg)) {
    m <- tb[tb$cfg == "mb", ]; t_ <- tb[tb$cfg == "tutte", ]
    cat(sprintf("P_F1: componente piu' grande, mediana (solo define_main_branches fedele): %.1f%% (soglia 95%%)\n", median(m$comp_grande_pct)))
    cat(sprintf("P_F2: variazione mediana del DBHcyl: %.1f punti (soglia 0.7)\n", median(abs(m$errF_nuovo - m$errF_v4))))
    cat(sprintf("P_F3: alberi con altezza QSM aumentata di almeno il 5%%: %d su %d (soglia 3 su 5)\n", sum(m$var_h_pct >= 5), nrow(m)))
    cat(sprintf("P_F4: tre funzioni fedeli: mediana |errore| DBHcyl %.1f%% (soglia 4.0%%)\n", median(abs(t_$errF_nuovo))))
  }
  invisible(tb)
}