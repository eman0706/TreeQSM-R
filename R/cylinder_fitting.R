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

#' Cylinder Fitting
#'
#' Porting fedele di cylinder_fitting.m (parte di cylinders.m v3.0.0,
#' Pasi Raumonen). Per segmenti con piu' di 6 layer, aggrega
#' adattivamente da 3 a 10 layer per regione (invece di fittare un
#' layer alla volta) cosi' ogni fit ha abbastanza punti per essere
#' affidabile; sceglie la lunghezza di regione migliore in base alla
#' copertura di superficie, poi avanza lungo il segmento aggiungendo
#' pesi diversi ai punti "di raccordo" tra una regione e la successiva.
#' Per segmenti piccoli (<=6 layer), fitta un solo cilindro su tutto il
#' segmento.
#'
#' @param P Point cloud (n x 3)
#' @param Points Indici dei punti nel segmento (vettore, concatenato per
#'   layer - vedi verticalcat())
#' @param Ind Matrice nl x 2: indice di inizio/fine, dentro Points, dei
#'   punti di ciascun layer
#' @param nl Numero di layer nel segmento
#' @param si Indice del segmento corrente (1 = tronco/radice)
#' @return Lista con cyl (i cilindri fittati: radius, length, start,
#'   axis, mad, SurfCov, conv, rel - vettori/matrici allineati per
#'   indice) e Reg (le coordinate dei punti usati per ciascun cilindro)
#' @export
cylinder_fitting <- function(P, Points, Ind, nl, si) {
  
  mrange <- function(a, b) if (b < a) integer(0) else a:b
  hgt <- function(Pts, base, axis) {
    if (is.null(dim(Pts))) Pts <- matrix(Pts, nrow = 1)
    as.vector(sweep(Pts, 2, base) %*% axis)
  }
  
  if (nl > 6L) {
    i0 <- 1L; i <- 4L
    t <- 0L
    Reg  <- vector("list", nl)
    cyls <- vector("list", 11)
    regs <- vector("list", 11)
    data <- matrix(0, 11, 4)
    cyl <- NULL
    field_names <- NULL
    CylTop <- NULL
    c0 <- list()   # persiste tra le iterazioni (fedele: le struct MATLAB
    # mantengono i campi non riassegnati esplicitamente,
    # a differenza di una list() R ricreata da zero)
    
    while (i0 < nl - 2L) {
      i0_prev <- i0  # per rilevare un'iterazione che non avanza, in caso di errore
      region_ok <- tryCatch({
        ## Fit at least three cylinders di lunghezza diversa
        bot <- Points[mrange(Ind[i0, 1], Ind[i0 + 1, 2])]
        Bot <- average(P[bot, , drop = FALSE])
        again <- TRUE
        j <- 0L
        
        while (i + j <= nl && j <= 10L && (j <= 2L || again)) {
          RegC <- Points[mrange(Ind[i0, 1], Ind[i + j, 2])]
          top  <- Points[mrange(Ind[i + j - 1, 1], Ind[i + j, 2])]
          Top  <- average(P[top, , drop = FALSE])
          Axis <- Top - Bot
          axis_len_j <- sqrt(sum(Axis^2))
          if (!is.finite(axis_len_j) || axis_len_j < 1e-6) {
            stop("asse quasi nullo (Top e Bot troppo vicini) - regione degenere")
          }
          c0$axis <- Axis / axis_len_j
          h <- hgt(P[RegC, , drop = FALSE], Bot, c0$axis)
          minh <- min(h)
          if (j == 0L) {
            Bot <- Bot + minh * c0$axis
            c0$start <- Bot
            h <- hgt(P[RegC, , drop = FALSE], Bot, c0$axis)
            minh <- min(h)
          }
          if (i + j >= nl) {
            ht <- hgt(matrix(Top, 1, 3), c0$start, c0$axis)
            Top <- Top + (max(h) - ht) * c0$axis
          }
          ht <- hgt(matrix(Top, 1, 3), c0$start, c0$axis)
          Sec <- h <= ht & h >= minh
          c0$length <- ht - minh
          reg <- RegC[Sec]
          Q0 <- P[reg, , drop = FALSE]
          
          ## Filtra i punti e stima il raggio
          if (nrow(Q0) > 20L) {
            if (Sys.getenv("CF_DEBUG2") == "1") {
              cat(sprintf("  [CF-DBG2] ramo grande: dim(c0$axis)=%s length=%d class=%s\n",
                          paste(dim(c0$axis), collapse="x"), length(c0$axis), class(c0$axis)))
            }
            scf <- surface_coverage_filtering(Q0, c0, 0.02, 20)
            reg <- reg[scf$Pass]
            Q0  <- Q0[scf$Pass, , drop = FALSE]
            c0  <- scf$c
          } else {
            c0$radius  <- 0.01
            c0$SurfCov <- 0.05
            c0$mad     <- 0.01
            c0$conv    <- TRUE
            c0$rel     <- TRUE
          }
          
          ## Fitta il cilindro
          if (nrow(Q0) > 9L) {
            if (i >= nl && t == 0L) {
              c <- least_squares_cylinder(Q0, c0)
            } else if (i >= nl && t > 0L) {
              # NOTA (vedi conversazione): in questo ramo l'originale MATLAB
              # applica il peso 1/3 al blocco n2 e 2/3 al blocco n1 - un
              # ordine che, isolato, sembra invertito rispetto alla logica
              # coerente degli altri due casi analoghi (dove il peso
              # maggiore va sempre al gruppo "sezione"/Q). Testato
              # empiricamente un fix che allinea i pesi in modo "coerente"
              # (1/3 a n1, 2/3 a n2, dando piu' peso a Q) - il risultato e'
              # stato una stabilita' nettamente PEGGIORE su tutto il batch
              # (735s invece di 168s, molti fallback a 1-2cm) rispetto al
              # comportamento fedele. Probabile che l'originale sia
              # calibrato insieme al resto delle soglie dell'algoritmo
              # (again/SurfCov/RL) in modi non ovvi - "correggere" un pezzo
              # isolato senza poter ricalibrare tutto il resto peggiora le
              # cose. NON ripetere questo esperimento senza un motivo
              # nuovo e diverso.
              hh <- hgt(Q0, CylTop, c0$axis)
              I  <- hh >= 0
              Q  <- Q0[I, , drop = FALSE]
              reg <- reg[I]
              n2 <- nrow(Q); n1 <- sum(!I)
              if (n2 > 9L && n1 > 5L) {
                Q0b <- rbind(Q0[!I, , drop = FALSE], Q)
                W   <- c(rep(1/3, n2), rep(2/3, n1))
                c <- least_squares_cylinder(Q0b, c0, W, Q)
              } else {
                c <- least_squares_cylinder(Q0, c0)
              }
            } else if (t == 0L) {
              top2 <- Points[mrange(Ind[i + j - 3, 1], Ind[i + j - 2, 2])]
              Top2 <- average(P[top2, , drop = FALSE])
              ht2  <- hgt(matrix(Top2, 1, 3), Bot, c0$axis)
              hh   <- hgt(Q0, Bot, c0$axis)
              I <- hh <= ht2
              Q <- Q0[I, , drop = FALSE]
              reg <- reg[I]
              n2 <- nrow(Q); n3 <- sum(!I)
              if (n2 > 9L && n3 > 5L) {
                Q0b <- rbind(Q, Q0[!I, , drop = FALSE])
                W   <- c(rep(2/3, n2), rep(1/3, n3))
                c <- least_squares_cylinder(Q0b, c0, W, Q)
              } else {
                c <- least_squares_cylinder(Q0, c0)
              }
            } else {
              top2 <- Points[mrange(Ind[i + j - 3, 1], Ind[i + j - 2, 2])]
              Top2 <- average(P[top2, , drop = FALSE])
              ht2  <- hgt(matrix(Top2, 1, 3), CylTop, c0$axis)
              hh   <- hgt(Q0, CylTop, c0$axis)
              I1 <- hh < 0
              I2 <- hh >= 0 & hh <= ht2
              I3 <- hh > ht2
              Q <- Q0[I2, , drop = FALSE]
              reg <- reg[I2]
              n1 <- sum(I1); n2 <- nrow(Q); n3 <- sum(I3)
              if (n2 > 9L) {
                Q0b <- rbind(Q0[I1, , drop = FALSE], Q, Q0[I3, , drop = FALSE])
                W   <- c(rep(1/4, n1), rep(2/4, n2), rep(1/4, n3))
                c <- least_squares_cylinder(Q0b, c0, W, Q)
              } else {
                c <- c0
                c$rel <- FALSE
              }
            }
            
            if (!isTRUE(c$conv)) { c <- c0; c$rel <- FALSE }
            if (!is.null(c$SurfCov) && is.finite(c$SurfCov) && c$SurfCov < 0.2) c$rel <- FALSE
            
          } else {
            c <- c0
            c$rel <- FALSE
          }
          
          RL <- if (is.finite(c$radius) && c$radius > 0) c$length / c$radius else 0
          data[j + 1, ] <- c(isTRUE(c$rel), isTRUE(c$conv), c$SurfCov %||% 0, RL)
          cyls[[j + 1]] <- c
          regs[[j + 1]] <- reg
          j <- j + 1L
          
          ## Se il cilindro e' gia' ragionevole, smetti di fittarne altri
          ## (ma fittane sempre almeno tre)
          if (again && isTRUE(c$rel) && isTRUE(c$conv) && RL > 2) {
            if (si == 1L && c$SurfCov > 0.7) {
              again <- FALSE
            } else if (si > 1L && c$SurfCov > 0.5) {
              again <- FALSE
            }
          }
        }
        
        ## Seleziona il migliore tra i cilindri fittati (max copertura)
        OKfit <- as.logical(data[seq_len(j), 1]) & as.logical(data[seq_len(j), 2]) &
          (data[seq_len(j), 4] > 1.5)
        J <- seq_len(j)
        t <- t + 1L
        if (any(OKfit)) J <- J[OKfit]
        score <- data[J, 3] - 0.01 * data[J, 4]
        best  <- J[which.max(score)]
        c <- cyls[[best]]
        
        ## Aggiorna gli indici dei layer per la prossima regione
        CylTop <- c$start + c$length * c$axis
        i0 <- i0 + 1L
        bot <- Points[mrange(Ind[i0, 1], Ind[i0 + 1, 2])]
        Bot <- average(P[bot, , drop = FALSE])
        hh  <- hgt(matrix(Bot, 1, 3), CylTop, c$axis)
        i00 <- i0
        while (i0 + 1 < nl && i0 < i00 + 5L && hh < -c$radius / 3) {
          i0 <- i0 + 1L
          bot <- Points[mrange(Ind[i0, 1], Ind[i0 + 1, 2])]
          Bot <- average(P[bot, , drop = FALSE])
          hh  <- hgt(matrix(Bot, 1, 3), CylTop, c$axis)
        }
        i <- min(i0 + 5L, nl)
        
        ## Se la sezione finale e' molto corta, allunga semplicemente
        ## il cilindro corrente invece di crearne un altro
        if (nl - i0 + 2 < 4) {
          reg <- Points[mrange(Ind[nl - 5, 1], Ind[nl, 2])]
          Q0  <- P[reg, , drop = FALSE]
          ht_end <- sum((c$start + c$length * c$axis) * c$axis)
          hh <- as.vector(Q0 %*% c$axis)
          maxh <- max(hh)
          if (maxh > ht_end) c$length <- c$length + (maxh - ht_end)
          i0 <- nl
        }
        Reg[[t]] <- regs[[best]]
        
        if (t == 1L) {
          cyl <- c
          cyl$start <- matrix(cyl$start, ncol = 3)
          cyl$axis  <- matrix(cyl$axis, ncol = 3)
          field_names <- names(cyl)
        } else {
          for (nm in field_names) {
            if (nm %in% c("start", "axis")) {
              cyl[[nm]] <- rbind(cyl[[nm]], c[[nm]])
            } else {
              cyl[[nm]] <- c(cyl[[nm]], c[[nm]])
            }
          }
        }
        
        CylTop <- c$start + c$length * c$axis
        TRUE  # regione completata con successo
      }, error = function(e) {
        cat(sprintf("  [CF-WARN] regione fallita (i0~%d): %s - uso quanto accumulato finora\n",
                    i0_prev, conditionMessage(e)))
        FALSE
      })
      if (!isTRUE(region_ok)) break
    }
    Reg <- Reg[seq_len(t)]
    
  } else {
    ## Segmento piccolo: un solo cilindro su tutto il segmento
    Q0 <- P[Points, , drop = FALSE]
    if (nrow(Q0) > 10L) {
      bot <- Points[mrange(Ind[1, 1], Ind[1, 2])]
      Bot <- average(P[bot, , drop = FALSE])
      top <- Points[mrange(Ind[nl, 1], Ind[nl, 2])]
      Top <- average(P[top, , drop = FALSE])
      Axis <- Top - Bot
      axis_len <- sqrt(sum(Axis^2))
      # FIX: se Top e Bot sono troppo vicini (segmento degenere, comune su
      # dati rumorosi vicino alla base), Axis/axis_len darebbe NaN,
      # propagandosi fino a orthonormal_vectors() e rompendo tutto a
      # valle. Scartiamo il segmento con grazia invece di propagare NaN.
      if (!is.finite(axis_len) || axis_len < 1e-6) {
        cyl <- NULL; t <- 0L; Reg <- list()
      } else {
        c0 <- list(axis = Axis / axis_len)
        hh <- as.vector(Q0 %*% c0$axis)
        c0$length <- max(hh) - min(hh)
        hpoint <- sum(Bot * c0$axis)
        c0$start <- Bot - (hpoint - min(hh)) * c0$axis
        
        scf <- surface_coverage_filtering(Q0, c0, 0.02, 20)
        Reg <- vector("list", 1)
        Reg[[1]] <- Points[scf$Pass]
        Q0 <- Q0[scf$Pass, , drop = FALSE]
        cyl <- least_squares_cylinder(Q0, scf$c)
        if (!isTRUE(cyl$conv) || !isTRUE(cyl$rel)) cyl <- scf$c
        # FIX: least_squares_cylinder() (e scf$c, il fallback) restituiscono
        # start/axis come VETTORI PIATTI lunghezza 3, non matrici 1x3 - a
        # differenza del ramo "segmento grande" dove l'accumulo con rbind()
        # garantisce sempre una vera matrice anche per un solo cilindro.
        # Senza questa normalizzazione, codice a valle (adjustments.R) che fa
        # cyl$start[1,] assumendo sempre una matrice fallisce con "numero di
        # dimensioni errato" ogni volta che un segmento produce un solo
        # cilindro tramite questo ramo.
        cyl$start <- matrix(cyl$start, 1, 3)
        cyl$axis  <- matrix(cyl$axis, 1, 3)
        t <- 1L
      }
    } else {
      cyl <- NULL
      t <- 0L
      Reg <- list()
    }
  }
  
  if (t > 0L && length(Reg) > 0L) {
    for (i in seq_len(t)) Reg[[i]] <- P[Reg[[i]], , drop = FALSE]
    Reg <- Reg[seq_len(t)]
  } else {
    Reg <- list()
  }
  
  # FIX: quando un segmento produce un solo cilindro, start/axis restano
  # vettori piatti (lunghezza 3) invece di matrici 1x3 - non passano mai
  # per rbind(), che li converte in matrice solo quando ci sono 2+
  # cilindri da accumulare. Qualunque codice a valle che fa cyl$start[1,]
  # (parent_cylinder, adjustments) aspettandosi una matrice va in crash
  # ("numero di dimensioni errato") su un vettore. Forziamo sempre la
  # forma a matrice Nx3 in uscita, indipendentemente da quale ramo
  # interno li ha prodotti.
  if (!is.null(cyl)) {
    if (is.null(dim(cyl$start))) cyl$start <- matrix(cyl$start, nrow = 1)
    if (is.null(dim(cyl$axis)))  cyl$axis  <- matrix(cyl$axis,  nrow = 1)
  }
  
  list(cyl = cyl, Reg = Reg)
}
