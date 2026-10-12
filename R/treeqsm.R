# ==============================================================================
# 1. CARICAMENTO MODULI (Metterli QUI)
# ==============================================================================
# I moduli vengono caricati da test_run.R prima di chiamare treeqsm()
# (source rimosso per evitare sovrascrittura delle patch)

#' Reconstructs quantitative structure tree models from point clouds
#' @param P Matrix (n x 3)
#' @param inputs List of parameters
#' @return A list of QSM models
#' @export
treeqsm <- function(P, inputs) {
  # 1. Validazione e Pulizia Input
  if (!is.matrix(P)) P <- as.matrix(P)
  P <- P[, 1:3]
  mode(P) <- "numeric"
  P <- P[!apply(P, 1, function(x) any(is.na(x))), ]
  
  # Setup Timing e Date
  start_time <- proc.time()
  Date <- matrix(0, 2, 6)
  Date[1, ] <- as.numeric(format(Sys.time(), c("%Y","%m","%d","%H","%M","%S")))
  
  # Parametri di default se mancanti
  if (is.null(inputs$BaseHeight)) inputs$BaseHeight <- 0.5
  if (is.null(inputs$TaperCor)) inputs$TaperCor <- 0
  if (is.null(inputs$Dist)) inputs$Dist <- FALSE
  
  # Estrazione parametri loop
  PatchDiam1 <- inputs$PatchDiam1
  PatchDiam2Min <- inputs$PatchDiam2Min
  PatchDiam2Max <- inputs$PatchDiam2Max
  BallRad1 <- inputs$BallRad1
  BallRad2 <- inputs$BallRad2
  
  if (inputs$disp > 0) {
    cat(sprintf("\n--- Avvio TreeQSM: %s | Albero: %d ---\n", inputs$name, inputs$tree))
  }# end of treeqsm function
  
  QSM <- list()
  nmodel <- 0
  
  # FIX: cover_sets() usa campionamento casuale non seedato
  # (sample.int/sample dentro generate_cover_uniform e
  # create_prioritized_permutation) - per costruzione dell'algoritmo
  # originale (produce covering diversi a runs diversi, come nel vero
  # MATLAB). Questo pero' rendeva impossibile distinguere l'effetto reale
  # di un fix dal semplice rumore campionario tra un lancio e l'altro
  # dello stesso identico albero. Seediamo in modo deterministico mA
  # albero-specifico: una firma dei dati (numero di punti + somma
  # coordinate) garantisce che lo STESSO identico albero produca sempre
  # lo stesso risultato run dopo run, mentre alberi diversi restano
  # comunque diversificati tra loro (nessun seed fisso identico per tutti).
  #
  # FIX (vedi conversazione): incluso anche inputs$model nel calcolo del
  # seed - nell'originale MATLAB, cover_sets.m usa randperm() (non
  # seedato) per l'ordine di scansione dei punti, che PRODUCE COVER
  # DIVERSI per run ripetuti con lo STESSO input (commento originale:
  # "produces different covers for the same inputs") - questo e'
  # esattamente il meccanismo che select_optimum.m/Nmodels=5 sono
  # progettati per sfruttare (media/selezione tra run genuinamente
  # diversi). Il nostro seed fisso per albero eliminava questa
  # variabilita', rendendo repliche con lo stesso input identiche tra
  # loro (nessun beneficio da Nmodels>1). Con model=1 (il valore usato
  # finora ovunque nella pipeline) il contributo e' ESATTAMENTE ZERO -
  # comportamento identico a prima, nessun rischio per la stabilita'
  # gia' raggiunta. Per model=2,3,4,5 (nuove repliche con lo stesso
  # input) il seed cambia, permettendo vera variabilita'.
  model_val <- if (!is.null(inputs$model)) as.integer(inputs$model) else 1L
  if (is.na(model_val) || length(model_val) == 0) model_val <- 1L
  seed_val <- as.integer(abs(round(sum(P[, 1]) * 1000 + sum(P[, 2]) * 997 +
                                     nrow(P) * 991 +
                                     (model_val - 1L) * 104729)) %% 2147483647L)
  if (seed_val == 0L) seed_val <- 1L
  set.seed(seed_val)
  
  for (h in seq_along(PatchDiam1)) {
    # Step 1: Generazione Cover Iniziale
    inputs$PatchDiam1 <- PatchDiam1[h]
    inputs$BallRad1   <- BallRad1[min(h, length(BallRad1))]
    cover1 <- cover_sets(P, inputs)
    
    # Step 2: Identificazione Albero (Tree Sets)
    ts1 <- tree_sets(P, cover1, inputs)
    # Fix robusto per Base/Forb
    Base1 <- if(is.null(ts1$Base)) which.min(as.matrix(cover1$center)[,3]) else ts1$Base
    
    # Step 3: Segmentazione Iniziale (per stima Relative Size)
    segment1 <- segments(ts1$cover, Base1, ts1$Forb)
    segment1 <- correct_segments(P, ts1$cover, segment1, inputs, 0, 1, 1)
    
    for (i in seq_along(PatchDiam2Max)) {
      inputs$PatchDiam2Max <- PatchDiam2Max[i]
      inputs$BallRad2      <- BallRad2[min(i, length(BallRad2))]
      
      for (j in seq_along(PatchDiam2Min)) {
        inputs$PatchDiam2Min <- PatchDiam2Min[j]
        
        # --- GENERAZIONE MODELLO FINALE ---
        # Calcolo Relative Size (RS) per adattare i patch ai rami
        RS <- relative_size(P, ts1$cover, segment1)
        
        # Seconda cover (Isotropica o quasi, per evitare errori di orientamento)
        cover2 <- cover_sets(P, inputs, RS)
        ts2 <- tree_sets(P, cover2, inputs, segment1)
        
        # --- FIX ROBUSTO BASE ---
        Base2 <- ts2$Base
        if (any(is.na(Base2)) || length(Base2) == 0 || length(cover2$center) == 0) {
          # Fallback: usa il centro con altezza minima
          if (length(cover2$center) > 0) {
            cen_coords <- P[cover2$center, , drop = FALSE]
            Base2 <- which.min(cen_coords[, 3])
          } else {
            Base2 <- 1
          }
        }
        Forb2 <- if(is.logical(ts2$Forb)) which(ts2$Forb) else ts2$Forb
        
        # Segmentazione Finale
        segment2 <- segments(ts2$cover, Base2, Forb2)
        segment2 <- correct_segments(P, ts2$cover, segment2, inputs, 1, 1, 0)
        
        if (length(segment2$segments) == 0) next
        
        # Fitting Cilindri
        cylinder <- cylinders(P, ts2$cover, segment2, inputs)
        
        if (!is.null(cylinder$radius) && length(cylinder$radius) > 0) {
          nmodel <- nmodel + 1
          
          # Calcoli Strutturali
          branch <- tryCatch(branches(cylinder), error = function(e) {
            cat("  WARN branches():", e$message, "\n")
            list(order=integer(length(cylinder$radius)),
                 parent=integer(length(cylinder$radius)),
                 volume=numeric(length(cylinder$radius)),
                 length=numeric(length(cylinder$radius)),
                 angle=numeric(length(cylinder$radius)),
                 height=numeric(length(cylinder$radius)),
                 azimuth=numeric(length(cylinder$radius)),
                 diameter=numeric(length(cylinder$radius)))
          })
          
          # Estrazione Tronco (per calcoli volume/altezza)
          # Il tronco è il segmento con il massimo numero di layer (coerente con cylinders.R)
          seg_nl_t    <- sapply(segment2$segments, length)
          # FIX (coerenza con branch==1, fedele): il tronco e' sempre il
          # segmento 1 (root della segmentazione) - non "quello con
          # piu' layer" (euristica residua di prima che modify_topology
          # garantisse che il tronco fisico sia catturato per intero dal
          # segmento 1, ora potenzialmente inconsistente con branch==1
          # usato ovunque altrove, es. in tree_data.R).
          trunk_si    <- 1L
          # trunk_idx = indici dei cover set del segmento tronco
          trunk_cs    <- unlist(segment2$segments[[trunk_si]])
          nb_cov      <- length(ts2$cover$ball)
          trunk_cs    <- trunk_cs[!is.na(trunk_cs) & trunk_cs >= 1L & trunk_cs <= nb_cov]
          trunk_idx   <- unique(unlist(ts2$cover$ball[trunk_cs]))
          trunk_idx   <- trunk_idx[!is.na(trunk_idx) & trunk_idx >= 1L & trunk_idx <= nrow(P)]
          if (length(trunk_idx) == 0L) trunk_idx <- seq_len(min(100L, nrow(P)))
          trunk_points <- P[trunk_idx, , drop = FALSE]
          
          # Calcoli Dati Albero (Volume, Altezza, etc.)
          tdata <- tree_data(cylinder, branch, trunk_points, inputs)
          
          # Distanze punto-modello (fedele: if inputs.Dist nell'originale)
          pmdis <- NULL
          if (isTRUE(inputs$Dist)) {
            pmdis <- tryCatch(
              point_model_distance(P, cylinder),
              error = function(e) {
                cat("  WARN point_model_distance():", conditionMessage(e), "\n")
                NULL
              }
            )
          }
          
          # Chiusura Modello
          Date[2, ] <- as.numeric(format(Sys.time(), c("%Y","%m","%d","%H","%M","%S")))
          
          QSM[[nmodel]] <- list(
            cylinder = cylinder,
            branch = branch,
            treedata = tdata$treedata,
            pmdistance = pmdis,
            rundata = list(
              inputs = inputs,
              time = (proc.time() - start_time)[3],
              date = Date,
              version = "0.1.1-R_Optimized"
            )
          )
          
          if (inputs$disp >= 1) {
            cat(sprintf(" Modello %d: %d cilindri | Vol: %.4f m3 | Alt: %.2f m\n", 
                        nmodel, length(cylinder$radius), 
                        tdata$treedata$Volume[1], tdata$treedata$Height[1]))
          }
        }
      }
    }
  }
  return(QSM)
}