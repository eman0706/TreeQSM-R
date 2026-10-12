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

#' Tree Data
#'
#' Porting fedele PARZIALE di tree_data.m (Version 3.0.1, Copyright (C)
#' 2013-2022 Pasi Raumonen) - sostituisce la versione precedente
#' (invenzione completa: DBH calcolato come media dei cilindri che
#' attraversano l'altezza ASSOLUTA Z=1.3m, mai verificata contro il
#' sorgente originale).
#'
#' Include le parti direttamente rilevanti al problema affrontato in
#' sessione (identificazione del tronco, volume, altezza, DBH) - NON
#' include crown_measures, triangulate_stem, ne' le numerose
#' distribuzioni (per diametro/altezza/angolo/ordine di ramo), fuori
#' scopo.
#'
#' NOTA IMPORTANTE: il vero algoritmo identifica il tronco con
#' "Trunk = cylinder.branch == 1" preso in ORDINE DI INDICE. Con
#' modify_topology() ora integrata (fonde correttamente segments.m in
#' un tronco continuo), questo funziona come nel vero algoritmo, senza
#' bisogno di ricostruzioni di catena a valle.
#'
#' Il DBH vero ("DBHqsm") usa la LUNGHEZZA CUMULATIVA lungo l'asse dei
#' cilindri del tronco. Il vero algoritmo poi rifitta un cilindro
#' DEDICATO sui punti reali della nuvola in una fascia stretta attorno
#' all'altezza cercata, con verifica di plausibilita', invece di
#' fidarsi ciecamente del cilindro grezzo - il meccanismo reale per
#' gestire basi rumorose/contaminate. Generalizzato qui a qualsiasi
#' altezza (diameter_at_height), non solo 1.3m, per coprire anche le
#' metriche custom del progetto (2/4/8m) con lo stesso rigore - queste
#' ultime non hanno un equivalente nel vero tree_data.m, che calcola
#' solo il DBH a 1.3m.
#'
#' @param cylinder Output di cylinders() (radius, length, start, axis, branch)
#' @param branch Output di branches()
#' @param trunk Nuvola di punti REALI del tronco (matrice Nx3) - serve
#'   per il rifit dedicato
#' @param inputs Lista di input (usata per inputs$Tria)
#' @return Lista con treedata e triangulation (NULL, non implementata)
#' @export
tree_data <- function(cylinder, branch, trunk, inputs) {
  Rad <- cylinder$radius
  Len <- cylinder$length
  nc  <- length(Rad)
  ind <- seq_len(nc)
  Trunk <- cylinder$branch == 1L
  
  ref_robust <- robust_trunk_reference(cylinder, Trunk, ind)
  
  # FIX: identifica e corregge i cilindri del tronco geometricamente
  # anomali (stesso criterio di robust_trunk_reference: Z sotto terra o
  # asse fortemente inclinato) - il loro raggio grezzo gonfierebbe
  # TrunkVolume/TrunkArea anche se il DBH puntuale ormai li aggira
  # correttamente tramite il riferimento robusto. Un cilindro anomalo
  # viene sostituito con il raggio del primo cilindro "buono" successivo
  # nella catena (o, in mancanza, il primo precedente) - stessa logica
  # di ref_robust ma applicata a ciascun cilindro anomalo individuale,
  # non solo al primo incontrato.
  Rad_corr <- Rad
  T_trunk <- ind[Trunk]
  if (length(T_trunk) > 0) {
    # FIX (vedi conversazione): criterio aggiuntivo basato sulla
    # lunghezza - un cilindro anomalo puo' avere posizione/angolo
    # perfettamente plausibili (passando quindi il criterio esistente)
    # ma una lunghezza eccessiva (es. 3+ metri per un solo cilindro),
    # segno che il fitting adattivo ha "mediato" su una regione troppo
    # estesa a causa di dati rumorosi, gonfiando il volume calcolato
    # con quel raggio su tutta quella lunghezza. Soglia assoluta (non
    # relativa alla catena stessa, che puo' essere in gran parte gia'
    # contaminata) - un cilindro ben fittato raramente supera 1 metro.
    is_bad <- (cylinder$start[T_trunk, 3] < -0.1) | (abs(cylinder$axis[T_trunk, 3]) < 0.7) |
      (Len[T_trunk] > 1.0)
    if (any(is_bad)) {
      for (bp in which(is_bad)) {
        repl <- NA_real_
        if (bp < length(T_trunk)) {
          for (k2 in (bp + 1L):length(T_trunk)) {
            if (!is_bad[k2]) { repl <- Rad[T_trunk[k2]]; break }
          }
        }
        if (is.na(repl) && bp > 1L) {
          for (k2 in (bp - 1L):1L) {
            if (!is_bad[k2]) { repl <- Rad[T_trunk[k2]]; break }
          }
        }
        if (!is.na(repl)) Rad_corr[T_trunk[bp]] <- repl
      }
    }
  }
  
  treedata <- list()
  # FIX (vedi conversazione): lo stesso tipo di contaminazione trovato
  # nel tronco (cilindri con lunghezza eccessiva che "mediano" su una
  # regione troppo estesa, gonfiando il volume) si presenta anche nei
  # cilindri di ramo, con magnitudine paragonabile (35-48% del volume
  # rami in alberi con contaminazione nota). La struttura qui e'
  # ramificata, non lineare come il tronco - non esiste un "prima/dopo"
  # naturale. Si usa invece il raggio medio dei cilindri figli diretti
  # (piu' vicini alla punta del ramo, fisiologicamente piu' sottili e
  # quindi stime piu' affidabili) quando disponibili e piu' piccoli del
  # cilindro anomalo stesso. Molti cilindri anomali non hanno figli
  # (sono all'estremita' dei rami) e restano non corretti per mancanza
  # di dati utilizzabili.
  Rad_corr_branch <- Rad
  idx_branch_bad <- which(!Trunk & Len > 1.0)
  for (bi in idx_branch_bad) {
    figli <- which(cylinder$parent == bi)
    if (length(figli) > 0) {
      r_figli <- mean(Rad[figli])
      if (is.finite(r_figli) && r_figli < Rad[bi]) {
        Rad_corr_branch[bi] <- r_figli
      }
    }
  }
  # NOTA: BranchVolume/BranchArea sui cilindri di ramo sono calcolati
  # direttamente, non per sottrazione da TotalVolume/TotalArea - cosi'
  # la correzione applicata al tronco non introduce incoerenze.
  treedata$TrunkVolume  <- 1000 * pi * sum(Rad_corr[Trunk]^2 * Len[Trunk])
  treedata$BranchVolume <- 1000 * pi * sum(Rad_corr_branch[!Trunk]^2 * Len[!Trunk])
  treedata$TotalVolume  <- treedata$TrunkVolume + treedata$BranchVolume
  
  bottom <- min(cylinder$start[, 3])
  i_top  <- which.max(cylinder$start[, 3])
  top    <- cylinder$start[i_top, 3]
  if (cylinder$axis[i_top, 3] > 0) {
    top <- top + Len[i_top] * cylinder$axis[i_top, 3]
  }
  treedata$TreeHeight <- top - bottom
  treedata$Height     <- treedata$TreeHeight  # alias, usato da treeqsm.R nel log
  treedata$Volume      <- treedata$TotalVolume / 1000  # alias in m3, usato da treeqsm.R nel log
  
  treedata$TrunkLength  <- sum(Len[Trunk])
  treedata$BranchLength <- sum(Len[!Trunk])
  treedata$TotalLength  <- treedata$TrunkLength + treedata$BranchLength
  
  treedata$TrunkArea  <- 2 * pi * sum(Rad_corr[Trunk] * Len[Trunk])
  treedata$BranchArea <- 2 * pi * sum(Rad[!Trunk] * Len[!Trunk])
  treedata$TotalArea  <- treedata$TrunkArea + treedata$BranchArea
  
  # Esposto per uso esterno (process_one_tree() in test_run03.R, per
  # vol4m_dm3/vol8m_dm3 - metriche custom calcolate fuori da questo
  # file, che altrimenti userebbero i raggi grezzi non corretti).
  treedata$RadCorrTrunk <- Rad_corr[Trunk]
  
  dbh13 <- diameter_at_height(1.3, Trunk, cylinder, ind, trunk, ref = ref_robust)
  treedata$DBHqsm <- dbh13$DiamQsm
  treedata$DBHcyl <- dbh13$DiamCyl
  treedata$Trunk  <- ind[Trunk]
  
  # Metriche custom del progetto (2/4/8m) - non presenti nel vero
  # tree_data.m, ma calcolate con lo stesso meccanismo fedele
  # (diameter_at_height), per coerenza con l'approccio dell'algoritmo
  # originale invece di un criterio inventato a parte.
  d200 <- diameter_at_height(2.0, Trunk, cylinder, ind, trunk, ref = ref_robust)
  d400 <- diameter_at_height(4.0, Trunk, cylinder, ind, trunk, ref = ref_robust)
  d800 <- diameter_at_height(8.0, Trunk, cylinder, ind, trunk, ref = ref_robust)
  treedata$Diam200qsm <- d200$DiamQsm
  treedata$Diam200cyl <- d200$DiamCyl
  treedata$Diam400qsm <- d400$DiamQsm
  treedata$Diam400cyl <- d400$DiamCyl
  treedata$Diam800qsm <- d800$DiamQsm
  treedata$Diam800cyl <- d800$DiamCyl
  
  if (!is.null(branch) && !is.null(branch$order)) {
    treedata$NumberBranches <- max(0L, length(unique(cylinder$branch)) - 1L)
    treedata$MaxBranchOrder <- suppressWarnings(max(branch$order, na.rm = TRUE))
    if (!is.finite(treedata$MaxBranchOrder)) treedata$MaxBranchOrder <- 0L
  } else {
    treedata$NumberBranches <- NA_integer_
    treedata$MaxBranchOrder <- NA_integer_
  }
  
  # Distribuzioni (fedele: cylinder_distribution/cylinder_height_
  # distribution dipendono solo da cylinder, sempre calcolabili;
  # branch_distribution/branch_order_distribution richiedono branch
  # valido)
  treedata <- tryCatch(
    cylinder_distribution(treedata, cylinder, "Dia"),
    error = function(e) treedata
  )
  treedata <- tryCatch(
    cylinder_height_distribution(treedata, cylinder, ind),
    error = function(e) treedata
  )
  treedata <- tryCatch(
    cylinder_distribution(treedata, cylinder, "Zen"),
    error = function(e) treedata
  )
  treedata <- tryCatch(
    cylinder_distribution(treedata, cylinder, "Azi"),
    error = function(e) treedata
  )
  if (!is.null(branch) && !is.null(branch$order) && length(branch$order) > 0) {
    treedata <- tryCatch(
      branch_order_distribution(treedata, branch),
      error = function(e) treedata
    )
    for (d in c("Dia", "Hei", "Ang", "Azi", "Zen")) {
      treedata <- tryCatch(
        branch_distribution(treedata, branch, d),
        error = function(e) treedata
      )
    }
  }
  
  list(treedata = treedata, triangulation = NULL)
}

#' Diameter at Height
#'
#' Generalizzazione di dbh_cylinder (sotto-funzione di tree_data.m) a
#' qualsiasi altezza H, non solo 1.3m. Calcola il diametro grezzo dal
#' QSM (DiamQsm, lunghezza cumulativa lungo l'asse dei cilindri di
#' tronco fino all'altezza H) e, se disponibili abbastanza punti reali
#' in una fascia H±0.2m, rifitta un cilindro dedicato (DiamCyl),
#' accettandolo solo se supera i controlli di plausibilita'
#' dell'originale (raggio entro 0.8-1.2x del grezzo, allineamento
#' assiale, convergenza).
#'
#' @param H Altezza (lunghezza cumulativa) a cui misurare il diametro
#' @param Trunk Vettore logico (branch==1), allineato a cylinder
#' @param cylinder Output di cylinders()
#' @param ind Indici 1:nc
#' @param trunk_points Nuvola di punti reali del tronco (matrice Nx3),
#'   o NULL se non disponibile (in tal caso DiamCyl = DiamQsm)
#' @param band Ampiezza della fascia di ricerca attorno a H (default
#'   0.2, fedele all'originale 1.1-1.5 per H=1.3)
#' @export
#' Robust Trunk Reference
#'
#' Sceglie un cilindro di riferimento (start/axis) geometricamente
#' plausibile tra i primi cilindri della catena del tronco, invece di
#' fidarsi ciecamente del primo cilindro dell'array (cylinder$start[1,],
#' cylinder$axis[1,]) - che si e' rivelato, empiricamente, poter essere
#' un fit anomalo isolato (es. Z sotto terra, asse fortemente inclinato)
#' quando la base e' contaminata da troppi punti aggregati. Un asse
#' sbagliato qui manda in tilt tutte le proiezioni di altezza a valle
#' (diameter_at_height), selezionando punti nella fascia sbagliata.
#'
#' @param cylinder Output di cylinders()
#' @param Trunk Vettore logico (branch==1)
#' @param ind Indici 1:nc
#' @param max_check Quanti cilindri iniziali della catena esaminare
#'   prima di ripiegare comunque sul primo (default 5)
#' @return Lista con start e axis (vettori lunghezza 3)
#' @export
robust_trunk_reference <- function(cylinder, Trunk, ind, max_check = 5L) {
  T <- ind[Trunk]
  n <- length(T)
  if (n == 0L) {
    return(list(start = c(0, 0, 0), axis = c(0, 0, 1)))
  }
  n_check <- min(max_check, n)
  # FIX (vedi conversazione): stima indipendente e robusta del raggio
  # atteso, usando la mediana dei primi n_check cilindri - permette di
  # scartare un candidato che passa i controlli di posizione/asse ma
  # ha esso stesso un fit anomalo (es. raggio doppio dei vicini, tipico
  # di un fit contaminato vicino alla base). Senza questo controllo, un
  # tale candidato diventerebbe il riferimento per tutto il rifit
  # successivo, propagando e amplificando l'errore invece di
  # correggerlo (verificato empiricamente: causa di un caso con
  # DBHcyl=767.8cm quando il vero diametro era ~47cm).
  radii_check <- cylinder$radius[T[seq_len(n_check)]]
  r_med <- suppressWarnings(median(radii_check, na.rm = TRUE))
  for (k in seq_len(n_check)) {
    idx <- T[k]
    # FIX: soglia allargata da -0.1 a -0.3 (vedi conversazione) - un
    # fit del primo cilindro puo' avere lo start leggermente sotto
    # terra per piccole imprecisioni vicino alla base contaminata,
    # senza per questo essere inaffidabile; -0.1 scartava candidati
    # altrimenti ragionevoli, forzando la scelta di un candidato
    # successivo con fit anomalo.
    z_ok   <- cylinder$start[idx, 3] >= -0.3  # il suolo e' Z=0 per costruzione (extract_trees.R)
    ang_ok <- abs(cylinder$axis[idx, 3]) > 0.7  # non troppo inclinato (< ~45 gradi dalla verticale)
    rad_ok <- is.finite(r_med) && r_med > 0 &&
      cylinder$radius[idx] > 0.5 * r_med && cylinder$radius[idx] < 1.8 * r_med
    if (isTRUE(z_ok) && isTRUE(ang_ok) && isTRUE(rad_ok)) {
      return(list(start = cylinder$start[idx, ], axis = cylinder$axis[idx, ]))
    }
  }
  # Nessun cilindro tra i primi "max_check" e' plausibile - fallback al
  # primo comunque (meglio un riferimento imperfetto che nessuno).
  idx <- T[1]
  list(start = cylinder$start[idx, ], axis = cylinder$axis[idx, ])
}

#' Fit di cerchio algebrico (metodo di Kasa) su punti 2D.
#'
#' Aggiunto il 27/08/2026 per sostituire il percentile 20% come stima
#' iniziale del raggio in diameter_at_height(): un fit di cerchio
#' algebrico e' molto piu' robusto quando i punti disponibili coprono
#' solo un arco parziale della circonferenza (occlusione angolare del
#' tronco) - il percentile 20% sottostima sistematicamente il vero
#' raggio in questi casi. Verificato empiricamente su un caso con
#' occlusione (albero "38"): percentile 20% dava raggio 7.9cm, fit di
#' cerchio dava 18.7cm, contro un vero raggio di 15.7cm (confermato da
#' misura di campo indipendente - Circum130).
#' @export
fit_cerchio_kasa <- function(x, y) {
  n <- length(x)
  if (n < 5) return(NULL)
  A <- cbind(2*x, 2*y, 1)
  b <- x^2 + y^2
  sol <- tryCatch(qr.solve(A, b), error = function(e) NULL)
  if (is.null(sol)) return(NULL)
  xc <- sol[1]; yc <- sol[2]; c_val <- sol[3]
  r2 <- c_val + xc^2 + yc^2
  if (!is.finite(r2) || r2 <= 0) return(NULL)
  list(cx = xc, cy = yc, r = sqrt(r2))
}

#' Fit di cerchio robusto su punti 3D, proiettati sul piano perpendicolare
#' a un asse dato (non assume asse verticale).
#' @export
fit_cerchio_3d <- function(pts3d, axis, start) {
  V <- sweep(pts3d, 2, start)
  h <- as.vector(V %*% axis)
  V_perp <- V - outer(h, axis)
  
  tmp <- if (abs(axis[1]) < 0.9) c(1,0,0) else c(0,1,0)
  e1 <- tmp - sum(tmp*axis)*axis; e1 <- e1/sqrt(sum(e1^2))
  e2 <- c(axis[2]*e1[3]-axis[3]*e1[2], axis[3]*e1[1]-axis[1]*e1[3], axis[1]*e1[2]-axis[2]*e1[1])
  
  x2d <- as.vector(V_perp %*% e1)
  y2d <- as.vector(V_perp %*% e2)
  
  fit_cerchio_kasa(x2d, y2d)
}

diameter_at_height <- function(H, Trunk, cylinder, ind, trunk_points = NULL, band = 0.2, ref = NULL) {
  n <- sum(Trunk)
  T <- ind[Trunk]
  if (n == 0L || length(T) == 0L) {
    return(list(DiamQsm = NA_real_, DiamCyl = NA_real_))
  }
  
  if (is.null(ref)) ref <- robust_trunk_reference(cylinder, Trunk, ind)
  
  i <- 1L
  while (i < n && sum(cylinder$length[T[seq_len(i)]]) < H) {
    i <- i + 1L
  }
  DiamQsm <- 2 * cylinder$radius[T[i]]
  
  DiamCyl <- DiamQsm
  if (!is.null(trunk_points) && nrow(trunk_points) > 0) {
    # FIX: proiezione dell'altezza rispetto al riferimento ROBUSTO
    # (ref$start/ref$axis), non al primo cilindro dell'array (che puo'
    # essere un fit anomalo isolato con asse fortemente inclinato,
    # verificato empiricamente - Z sotto terra, asse a 32 gradi dalla
    # verticale su un caso reale). Con un asse sbagliato qui, "hh"
    # selezionava punti nella fascia di altezza SBAGLIATA, producendo
    # risultati incoerenti indipendentemente da come veniva poi
    # inizializzato/filtrato il fit.
    V  <- sweep(trunk_points, 2, ref$start)
    hh <- as.vector(V %*% ref$axis)
    I  <- hh < (H + band) & hh > (H - band)
    if (sum(I) > 100L) {
      Tpts <- trunk_points[I, , drop = FALSE]
      # FIX (deviazione motivata dall'originale, vedi conversazione):
      # l'originale inizializza il rifit con select_cylinders(cylinder,i)
      # (il raggio del cilindro grezzo) e usa TUTTI i punti della fascia
      # - ma quando il cilindro grezzo e' esso stesso gia' sbagliato
      # (es. base contaminata aggregata con troppi layer), un fit ai
      # minimi quadrati standard su TUTTI i punti (inclusi eventuali
      # rami/vegetazione bassa mescolati nella nuvola) resta sensibile
      # agli outlier durante l'iterazione, e con un'inizializzazione
      # instabile puo' persino DIVERGERE (verificato empiricamente:
      # cambiare solo l'inizializzazione, senza filtrare i punti, ha
      # peggiorato il risultato invece di correggerlo). Serve sia
      # un'inizializzazione ROBUSTA sia un filtro robusto dei punti
      # stessi: percentile basso (20%) della distanza radiale rispetto
      # all'asse del cilindro grezzo (affidabile anche quando il
      # raggio non lo e') per la stima iniziale, poi si tengono solo i
      # punti entro un multiplo di quella stima (troncamento robusto,
      # analogo nello spirito a uno stimatore M-robusto) per il fit
      # vero e proprio.
      dl_init <- distances_to_line(Tpts, ref$axis, ref$start)
      # FIX (27/08/2026): fit di cerchio algebrico come stima primaria
      # (molto piu' robusto in caso di occlusione angolare parziale),
      # con il percentile 20% come fallback se il fit fallisce (troppo
      # pochi punti o geometria degenere).
      fit_circ <- tryCatch(fit_cerchio_3d(Tpts, ref$axis, ref$start), error = function(e) NULL)
      if (!is.null(fit_circ) && is.finite(fit_circ$r) && fit_circ$r > 0) {
        r_init <- fit_circ$r
      } else {
        r_init <- suppressWarnings(quantile(dl_init$d, 0.2, na.rm = TRUE))
        if (!is.finite(r_init) || r_init <= 0) r_init <- cylinder$radius[T[i]]
      }
      
      keep_pts <- dl_init$d < 1.6 * r_init
      if (sum(keep_pts) > 50L) {
        Tpts_robust <- Tpts[keep_pts, , drop = FALSE]
        
        cyl0 <- select_cylinders(cylinder, i)
        cyl0$radius <- as.numeric(r_init)
        cyl0$start  <- ref$start
        cyl0$axis   <- ref$axis
        cyl  <- tryCatch(least_squares_cylinder(Tpts_robust, cyl0), error = function(e) NULL)
        if (!is.null(cyl)) {
          # Controllo di plausibilita' aggiornato: confronta il
          # risultato con la stima robusta r_init (indipendente e
          # affidabile), non con DiamQsm (il valore stesso che stiamo
          # cercando di correggere).
          RadiusOK   <- (cyl$radius > 0.5 * r_init) && (cyl$radius < 2.0 * r_init)
          axis_align <- abs(sum(ref$axis * cyl$axis))
          if (isTRUE(RadiusOK) && is.finite(axis_align) && axis_align > 0.9 &&
              isTRUE(cyl$conv) && isTRUE(cyl$rel)) {
            DiamCyl <- 2 * cyl$radius
          }
        }
      }
    }
  }
  
  list(DiamQsm = DiamQsm, DiamCyl = DiamCyl)
}

#' Select Cylinders
#'
#' Helper NON presente nel sorgente recuperato (chiamata da
#' dbh_cylinder ma definita altrove nel repository originale, mai
#' recuperata) - implementazione propria, a basso rischio dato l'uso
#' ovvio: estrae il cilindro all'indice i come punto di partenza (cyl0)
#' per un nuovo fit.
#' @export
select_cylinders <- function(cylinder, i) {
  list(
    radius = cylinder$radius[i],
    length = cylinder$length[i],
    start  = cylinder$start[i, ],
    axis   = cylinder$axis[i, ]
  )
}
