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

#' Fit cylinders to tree segments
#' Versione 3.1 - Integra il porting fedele di cylinder_fitting.m
#' (fitting adattivo a piu' regioni, 3-10 layer, invece di un layer alla
#' volta), parent_cylinder.m (vera geometria di attraversamento
#' asse-cilindro), e adjustments.m (correzioni finali di rastremazione,
#' verificata riga per riga come completa e fedele - vedi conversazione).
#' La costruzione di Points/IndPoints segue fedelmente cylinders.m
#' (verticalcat + BallSize/IndPoints).
#'
#' "branch" resta assegnato come si (fedele al vero comportamento
#' MATLAB nei nostri dati - VERIFICATO EMPIRICAMENTE, non solo per
#' ipotesi: cylinder.extension(PC)==c non si e' mai verificato in un
#' batch completo di 82 alberi con diagnostica attiva, quindi il ramo
#' "eredita branch dal genitore" dell'originale non si attiverebbe mai
#' comunque nel nostro contesto). La vera catena del tronco va
#' ricostruita a valle con identify_stem_chain().
#' Filter Layer Outliers
#'
#' Filtro spaziale robusto, layer per layer, contro contaminazione da
#' vegetazione bassa/sottobosco erroneamente inclusa nella nuvola di
#' punti durante l'estrazione LiDAR (segmentazione individuale-albero
#' imperfetta, es. li2012 - verificato empiricamente: il primo layer
#' del tronco di un albero specifico si estendeva su 5-7 metri in X/Y,
#' quando un vero tronco dovrebbe stare in 20-40cm).
#'
#' NOTA DI PROGETTAZIONE (dopo DUE tentativi precedenti):
#' 1) Un filtro basato sulla mediana/MAD del solo layer corrente non
#'    basta quando la contaminazione e' massiccia (i punti spuri sono
#'    la MAGGIORANZA del layer) - la mediana stessa risulta gia'
#'    "avvelenata".
#' 2) Applicare il filtro SEMPRE, indiscriminatamente su tutti gli
#'    alberi, e' pericoloso: su alberi con un tronco genuinamente
#'    grande, taglia via parte del vero tronco insieme alla
#'    contaminazione (verificato: la stabilita' del batch e' crollata
#'    da 74 a 37 alberi OK, con nuovi valori assurdi comparsi).
#' Soluzione: il filtro si attiva SOLO se c'e' evidenza inequivocabile
#' di contaminazione massiccia nel primo layer (una frazione
#' significativa di punti a distanza clamorosamente implausibile dal
#' centro atteso, es. oltre 2 metri - nessun vero tronco, per quanto
#' grande, li raggiunge mai) - altrimenti il layer resta intatto,
#' senza rischio per alberi normali o genuinamente grandi.
#'
#' @param P Point cloud completa (Nx3)
#' @param Points Vettore di indici punto, concatenato per layer
#' @param IndPoints Matrice nl x 2 (range di Points per ciascun layer)
#' @param nl Numero di layer
#' @param init_center Centro iniziale per il primo layer (default
#'   c(0,0), il centroide dell'albero)
#' @param radius_max Raggio massimo assoluto plausibile, usato SOLO se
#'   il filtro si attiva (default 1.0m - permissivo, per non tagliare
#'   alberi genuinamente grandi anche quando il filtro interviene)
#' @param contamination_check_radius Soglia di distanza "clamorosamente
#'   implausibile" per il controllo preliminare (default 2.0m)
#' @param contamination_frac_threshold Frazione minima di punti oltre
#'   contamination_check_radius perche' il filtro si attivi (default
#'   0.15, cioe' 15%)
#' @return Vettore logico (stessa lunghezza di Points), TRUE = tenuto
#' @export
filter_layer_outliers <- function(P, Points, IndPoints, nl,
                                  init_center = c(0, 0), radius_max = 1.0,
                                  contamination_check_radius = 2.0,
                                  contamination_frac_threshold = 0.15) {
  keep_mask <- rep(TRUE, length(Points))
  
  if (IndPoints[1, 2] < IndPoints[1, 1]) return(keep_mask)
  idx0 <- IndPoints[1, 1]:IndPoints[1, 2]
  if (length(idx0) < 10) return(keep_mask)
  pts0 <- P[Points[idx0], 1:2, drop = FALSE]
  dists0 <- sqrt(rowSums(sweep(pts0, 2, init_center)^2))
  frac_extreme <- mean(dists0 > contamination_check_radius)
  
  if (frac_extreme < contamination_frac_threshold) {
    return(keep_mask)
  }
  
  center <- init_center
  for (j in seq_len(nl)) {
    if (IndPoints[j, 2] < IndPoints[j, 1]) next
    idx_range <- IndPoints[j, 1]:IndPoints[j, 2]
    if (length(idx_range) < 5) next
    pt_idx <- Points[idx_range]
    pts <- P[pt_idx, 1:2, drop = FALSE]
    
    dists <- sqrt(rowSums(sweep(pts, 2, center)^2))
    med_d <- median(dists)
    mad_d <- median(abs(dists - med_d))
    thresh <- min(max(0.15, med_d + 5 * mad_d), radius_max)
    bad <- dists > thresh
    
    kept_idx <- idx_range[!bad]
    if (length(bad) > 0 && any(bad)) {
      keep_mask[idx_range[bad]] <- FALSE
    }
    if (length(kept_idx) >= 3) {
      pts_kept <- P[Points[kept_idx], 1:2, drop = FALSE]
      center <- apply(pts_kept, 2, median)
    }
  }
  keep_mask
}

cylinders <- function(P, cover, segment, inputs) {
  
  mrange <- function(a, b) if (b < a) integer(0) else a:b
  
  Segs <- segment$segments
  SPar <- segment$ParentSegment
  SChi <- segment$ChildSegment
  NumOfSeg <- length(Segs)
  
  if (NumOfSeg == 0) return(NULL)
  
  if (is.null(dim(SPar))) {
    # SPar e' un vettore semplice (formato fedele di correct_segments():
    # solo i genitori, senza la colonna dell'indice di layer - che
    # cylinders() non usa comunque mai). matrix(SPar, ncol=2) lo
    # "piegava" in modo scorretto (dividendo il vettore a meta') invece
    # di trattarlo come colonna dei genitori.
    SPar <- cbind(as.integer(SPar), 0L)
  }
  if (ncol(SPar) < 2)     SPar <- cbind(SPar, 0L)
  SPar[is.na(SPar)] <- 0L
  
  max_cyl <- max(as.integer(NumOfSeg) * 50L, 4000L)
  cylinder <- list(
    radius      = numeric(max_cyl),
    length      = numeric(max_cyl),
    start       = matrix(0, max_cyl, 3),
    axis        = matrix(0, max_cyl, 3),
    parent      = integer(max_cyl),
    extension   = integer(max_cyl),
    SurfCov     = numeric(max_cyl),
    branch      = integer(max_cyl),
    BranchOrder = integer(max_cyl),
    UnmodRadius = numeric(max_cyl),
    mad         = numeric(max_cyl)
  )
  c_idx <- 1L
  
  # BFS dei segmenti (fedele: bases = segmenti senza genitore)
  bases <- which(SPar[, 1] == 0L)
  if (length(bases) == 0) bases <- 1L
  
  SegmentIndex <- integer(0)
  queue   <- bases
  visited <- logical(NumOfSeg)
  visited[bases] <- TRUE
  while (length(queue) > 0) {
    curr  <- queue[1]; queue <- queue[-1]
    SegmentIndex <- c(SegmentIndex, curr)
    children <- unlist(SChi[curr])
    children <- children[!is.na(children) & children > 0L & children <= NumOfSeg]
    for (ch in children) {
      if (!visited[ch]) { visited[ch] <- TRUE; queue <- c(queue, ch) }
    }
  }
  if (length(SegmentIndex) == 0) return(NULL)
  
  root <- SegmentIndex[1]
  cat(sprintf("  [CYL] root_seg=%d (nl=%d)\n", root, length(Segs[[root]])))
  
  brad2 <- inputs$BallRad2 %||% inputs$BallRad1 %||% 0.40
  
  CiS <- vector("list", NumOfSeg)
  last_cyl_of_seg <- integer(NumOfSeg)
  
  for (k in seq_along(SegmentIndex)) {
    si <- SegmentIndex[k]
    if (si <= 0) next
    tryCatch({
      
      Seg <- Segs[[si]]
      nl  <- length(Seg)
      vc  <- verticalcat(Seg)
      Sets <- vc$Sets
      IndSets <- vc$IndSets
      ns <- length(Sets)
      
      Sets_valid <- Sets >= 1L & Sets <= length(cover$ball) & !is.na(Sets)
      BallSize <- integer(length(Sets))
      BallSize[Sets_valid] <- sapply(Sets[Sets_valid], function(s) length(cover$ball[[s]]))
      
      Points <- unlist(lapply(seq_along(Sets), function(ii) {
        if (Sets_valid[ii]) cover$ball[[Sets[ii]]] else integer(0)
      }))
      np <- length(Points)
      
      IndPoints <- matrix(1L, nl, 2)
      for (j in seq_len(nl)) {
        IndPoints[j, 2] <- sum(BallSize[mrange(IndSets[j, 1], IndSets[j, 2])])
      }
      IndPoints[, 2] <- cumsum(IndPoints[, 2])
      if (nl > 1) {
        IndPoints[2:nl, 1] <- IndPoints[2:nl, 1] + IndPoints[1:(nl - 1), 2]
      }
      
      # FIX (riattivazione sperimentale, vedi conversazione): condizionato
      # a FILTER_ENABLE per testarlo isolatamente prima di renderlo
      # permanente, dato lo storico complicato di questo filtro (3
      # tentativi precedenti, l'ultimo con soglia condizionale 15% che
      # scattava comunque su alberi sani). Ora relative_size e' fedele -
      # la segmentazione risultante potrebbe comportarsi diversamente.
      if (Sys.getenv("FILTER_ENABLE") == "1") {
        if (si == root) {
          keep_mask <- filter_layer_outliers(P, Points, IndPoints, nl, init_center = c(0, 0))
        } else {
          first_valid <- which(IndPoints[, 2] >= IndPoints[, 1])[1]
          init_c <- c(0, 0)
          if (!is.na(first_valid) && IndPoints[first_valid, 2] >= IndPoints[first_valid, 1]) {
            idx0 <- IndPoints[first_valid, 1]:IndPoints[first_valid, 2]
            if (length(idx0) >= 3) {
              init_c <- apply(P[Points[idx0], 1:2, drop = FALSE], 2, median)
            }
          }
          keep_mask <- filter_layer_outliers(P, Points, IndPoints, nl, init_center = init_c)
        }
        if (any(!keep_mask)) {
          new_counts <- integer(nl)
          new_points_list <- vector("list", nl)
          for (j in seq_len(nl)) {
            if (IndPoints[j, 2] < IndPoints[j, 1]) {
              new_points_list[[j]] <- integer(0)
              next
            }
            idx_range <- IndPoints[j, 1]:IndPoints[j, 2]
            kept <- idx_range[keep_mask[idx_range]]
            new_points_list[[j]] <- Points[kept]
            new_counts[j] <- length(kept)
          }
          Points <- unlist(new_points_list)
          np <- length(Points)
          IndPoints <- matrix(1L, nl, 2)
          IndPoints[, 2] <- cumsum(new_counts)
          if (nl > 1) {
            IndPoints[2:nl, 1] <- IndPoints[2:nl, 1] + IndPoints[1:(nl - 1), 2]
          }
        }
      }
      
      Base <- Seg[[1]]
      nb_base <- IndPoints[1, 2]
      
      # Ricostruisci solo segmenti abbastanza grandi (fedele)
      if (nl > 1 && np > nb_base && ns > 2 && np > 20 && length(Base) > 0) {
        
        fit_res <- tryCatch(
          cylinder_fitting(P, Points, IndPoints, nl, si),
          error = function(e) { cat(sprintf("  [CYL-ERR] si=%d cylinder_fitting: %s\n", si, conditionMessage(e))); NULL }
        )
        
        if (!is.null(fit_res) && !is.null(fit_res$cyl)) {
          cyl <- fit_res$cyl
          nc  <- length(cyl$radius)
          
          # FIX: rete di sicurezza esplicita PRIMA di adjustments() (ora
          # integrata e verificata come fedele, vedi conversazione).
          # Protezione aggiuntiva: un fit che supera i controlli interni di
          # cylinder_fitting (SurfCov, convergenza) ma cade su punti
          # isolati/rumorosi puo' comunque produrre un raggio fisicamente
          # assurdo (visto: 3-7 metri di raggio su un tronco vero) PRIMA
          # che adjustments() abbia la possibilita' di correggerlo.
          # Scartiamo le singole righe fuori da un limite fisico
          # ragionevole, non l'intero segmento.
          if (nc > 0) {
            radius_cap <- min(2.0, brad2 * 5)
            keep <- is.finite(cyl$radius) & cyl$radius > 0 & cyl$radius <= radius_cap
            if (!all(keep)) {
              cyl$radius  <- cyl$radius[keep]
              cyl$length  <- cyl$length[keep]
              cyl$start   <- cyl$start[keep, , drop = FALSE]
              cyl$axis    <- cyl$axis[keep, , drop = FALSE]
              cyl$mad     <- cyl$mad[keep]
              cyl$SurfCov <- cyl$SurfCov[keep]
              if (!is.null(cyl$conv)) cyl$conv <- cyl$conv[keep]
              if (!is.null(cyl$rel))  cyl$rel  <- cyl$rel[keep]
              nc <- length(cyl$radius)
            }
          }
          
          PC <- integer(0)
          added <- FALSE
          if (nc > 0 && si != root) {
            par_seg <- SPar[si, 1]
            if (par_seg > 0L && last_cyl_of_seg[par_seg] > 0L) {
              pcyl <- tryCatch(
                parent_cylinder(SPar, SChi, CiS, cylinder, cyl, si),
                error = function(e) NULL
              )
              if (!is.null(pcyl)) {
                PC <- pcyl$PC
                cyl <- pcyl$cyl
                added <- pcyl$added
                nc <- length(cyl$radius)
              }
            }
          }
          
          if (nc > 0) {
            cyl$radius0 <- cyl$radius
            
            # Porting fedele (parziale) di adjustments.m: MaxR/MinR dal
            # genitore, correzione di rastremazione, continuita' tra
            # cilindri. Sostituisce/integra il tetto fisico improvvisato
            # sopra con il vero meccanismo dell'algoritmo originale.
            parcyl <- if (length(PC) == 1 && PC > 0) {
              list(radius = cylinder$radius[PC], length = cylinder$length[PC],
                   start  = matrix(cylinder$start[PC, ], 1, 3),
                   axis   = matrix(cylinder$axis[PC, ], 1, 3))
            } else {
              list(radius = numeric(0), length = numeric(0),
                   start = matrix(0, 0, 3), axis = matrix(0, 0, 3))
            }
            cyl <- tryCatch(
              adjustments(cyl, parcyl, inputs, fit_res$Reg),
              error = function(e) { cat(sprintf("  [ADJ-ERR] si=%d: %s\n", si, conditionMessage(e))); cyl }
            )
            nc <- length(cyl$radius)
            
            Accept <- nc > 0 && min(cyl$radius[seq_len(nc)]) > 0
            if (Accept) {
              branch_id <- si   # fedele: vedi nota in cima al file
              
              start_idx <- c_idx
              for (kk in seq_len(nc)) {
                if (c_idx > max_cyl) break
                cylinder$radius[c_idx]      <- cyl$radius[kk]
                cylinder$length[c_idx]      <- cyl$length[kk]
                cylinder$start[c_idx, ]     <- cyl$start[kk, ]
                cylinder$axis[c_idx, ]      <- cyl$axis[kk, ]
                cylinder$UnmodRadius[c_idx] <- cyl$radius0[kk]
                cylinder$SurfCov[c_idx]     <- cyl$SurfCov[kk]
                cylinder$mad[c_idx]         <- cyl$mad[kk]
                cylinder$branch[c_idx]      <- branch_id
                
                if (kk == 1L) {
                  if (length(PC) == 1 && PC > 0) {
                    cylinder$parent[c_idx] <- PC
                    if (cylinder$extension[PC] == 0L) cylinder$extension[PC] <- c_idx
                  }
                  # altrimenti parent resta 0 (vera radice)
                } else {
                  cylinder$parent[c_idx]         <- c_idx - 1L
                  cylinder$extension[c_idx - 1L] <- c_idx
                }
                c_idx <- c_idx + 1L
              }
              end_idx <- c_idx - 1L
              if (end_idx >= start_idx) {
                CiS[[si]] <- start_idx:end_idx
                last_cyl_of_seg[si] <- end_idx
              }
            }
          }
        }
      }
    }, error = function(e) {
      cat(sprintf("  [CYL-LOOP-ERR] si=%d: %s\n", si, conditionMessage(e)))
      NULL
    })
    if (c_idx > max_cyl) break
  }
  
  total <- c_idx - 1L
  if (total == 0L) return(NULL)
  valid <- seq_len(total)
  
  res <- list(
    radius      = cylinder$radius[valid],
    length      = cylinder$length[valid],
    start       = cylinder$start[valid, , drop = FALSE],
    axis        = cylinder$axis[valid, , drop = FALSE],
    parent      = cylinder$parent[valid],
    extension   = cylinder$extension[valid],
    SurfCov     = cylinder$SurfCov[valid],
    branch      = cylinder$branch[valid],
    BranchOrder = integer(total),
    UnModRadius = cylinder$UnmodRadius[valid],
    mad         = cylinder$mad[valid],
    added       = logical(total),
    PositionInBranch = integer(total)
  )
  
  for (ii in valid) {
    p <- res$parent[ii]
    if (!is.na(p) && p > 0L && p <= total) {
      ext_p <- res$extension[p]
      res$BranchOrder[ii] <- res$BranchOrder[p] +
        if (!is.na(ext_p) && ext_p == ii) 0L else 1L
    }
  }
  for (b in unique(res$branch)) {
    idx_b <- which(res$branch == b)
    res$PositionInBranch[idx_b] <- seq_along(idx_b)
  }
  
  n_root_cyl <- sum(res$branch == root)
  cat(sprintf("  [CYL] %d cilindri | segmento radice=%d | branch unici=%d\n",
              total, n_root_cyl, length(unique(res$branch))))
  
  # Correzione growth volume (fedele: dentro cylinders.m stessa, non
  # nell'orchestratore principale) - opzionale, disattivata di default
  if (isTRUE(inputs$GrowthVolCor) && total > 0) {
    res <- tryCatch(
      growth_volume_correction(res, inputs),
      error = function(e) {
        cat("  WARN growth_volume_correction():", conditionMessage(e), "\n")
        res
      }
    )
  }
  
  return(res)
}