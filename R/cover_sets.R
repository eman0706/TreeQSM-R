#' Creates cover sets (surface patches) and neighbor relations for a point cloud
#'
#' @param P Point cloud matrix (n x 3)
#' @param inputs Input structure with fields: PatchDiam1, PatchDiam2Min, PatchDiam2Max, BallRad1, BallRad2, nmin1, nmin2
#' @param RelSize Optional relative sizes for variable-size cover
#' @return List with components: ball, center, neighbor
#' @export
cover_sets <- function(P, inputs, RelSize = NULL) {
  np <- nrow(P)
  
  # Make sure P is numeric matrix
  P <- as.matrix(P)
  mode(P) <- "numeric"
  
  # Determine if using uniform or variable size
  use_uniform <- is.null(RelSize)
  
  if (use_uniform) {
    # Uniform cover sets
    BallRad <- as.double(inputs$BallRad1)
    PatchDiamMax <- as.double(inputs$PatchDiam1)
    nmin <- as.double(inputs$nmin1)
    
    # Partition point cloud into cubes
    partition_result <- cubical_partition(P, BallRad)
    Partition <- partition_result$Partition
    CC <- partition_result$CubeCoord
    
    # Generate cover sets
    result <- generate_cover_uniform(P, Partition, CC, BallRad, PatchDiamMax, nmin)
  } else {
    # Variable size cover sets
    BallRad <- as.double(inputs$BallRad2)
    PatchDiamMin <- as.double(inputs$PatchDiam2Min)
    PatchDiamMax <- as.double(inputs$PatchDiam2Max)
    nmin <- as.double(inputs$nmin2)
    
    MRS <- PatchDiamMin / PatchDiamMax
    r <- 1.5 * (min(RelSize) / 256 * (1 - MRS) + MRS) * BallRad + 1e-5
    NE <- 1 + ceiling(BallRad / r)
    if (NE > 4) {
      r <- PatchDiamMax / 4
      NE <- 1 + ceiling(BallRad / r)
    }
    
    partition_result <- cubical_partition_with_cubes(P, r, NE)
    Partition <- partition_result$Partition
    CC <- partition_result$CubeCoord
    Cubes <- partition_result$Cubes
    
    # Mark points with no size
    RelSize[RelSize == 0] <- NA
    
    result <- generate_cover_variable(P, Partition, CC, Cubes, RelSize, 
                                      BallRad, PatchDiamMin, PatchDiamMax, 
                                      nmin, MRS, r)
  }
  
  return(result)
}

#' Generate uniform cover sets
generate_cover_uniform <- function(P, Partition, CC, BallRad, PatchDiamMax, nmin) {
  np <- nrow(P)
  Radius <- BallRad^2
  MaxDist <- PatchDiamMax^2
  
  # Initialize
  Ball <- list()
  Cen <- vector("integer", 0)
  NotExa <- rep(TRUE, np)
  Dist <- rep(1e8, np)
  BoP <- integer(np)
  nb <- 0
  
  # Random permutation
  RandPerm <- sample.int(np)
  
  # Get cube dimensions
  dims <- dim(Partition)
  
  for (i in 1:np) {
    Q <- RandPerm[i]
    if (NotExa[Q]) {
      # Get neighboring cubes
      cx <- CC[Q, 1]
      cy <- CC[Q, 2]
      cz <- CC[Q, 3]
      
      points <- c()
      for (dx in -1:1) {
        for (dy in -1:1) {
          for (dz in -1:1) {
            nx <- cx + dx
            ny <- cy + dy
            nz <- cz + dz
            if (nx >= 1 && nx <= dims[1] && 
                ny >= 1 && ny <= dims[2] && 
                nz >= 1 && nz <= dims[3]) {
              cell <- Partition[nx, ny, nz]
              if (!is.null(cell) && length(cell[[1]]) > 0) {
                points <- c(points, cell[[1]])
              }
            }
          }
        }
      }
      
      if (length(points) == 0) next
      
      # Make sure points is a valid integer vector
      if (is.list(points)) {
        points <- unlist(points)
      }
      points <- as.integer(points)
      points <- points[!is.na(points)]
      points <- points[points > 0 & points <= np]
      points <- unique(points)
      
      if (length(points) == 0) next
      
      # Compute distances
      V <- P[points, , drop = FALSE] - matrix(P[Q, ], nrow = length(points), ncol = 3, byrow = TRUE)
      dist <- rowSums(V^2)
      
      # Select points inside ball
      Inside <- dist < Radius
      if (sum(Inside) >= nmin) {
        ball <- points[Inside]
        d <- dist[Inside]
        core <- d < MaxDist
        NotExa[ball[core]] <- FALSE
        
        # Define new ball
        nb <- nb + 1
        Ball[[nb]] <- ball
        Cen[nb] <- Q
        
        # Update distances
        D <- Dist[ball]
        closer <- d < D
        ball <- ball[closer]
        Dist[ball] <- d[closer]
        BoP[ball] <- nb
      }
    }
  }
  
  # Trim vectors
  if (nb > 0) {
    Ball <- Ball[1:nb]
    Cen <- Cen[1:nb]
  }
  
  # Build cover sets
  cover <- build_cover_structure(P, Ball, Cen, nb, BoP)
  
  return(cover)
}

#' Generate variable-size cover sets
generate_cover_variable <- function(P, Partition, CC, Cubes, RelSize, 
                                    BallRad, PatchDiamMin, PatchDiamMax, 
                                    nmin, MRS, r) {
  np <- nrow(P)
  e <- BallRad - PatchDiamMax
  
  # Initialize
  Ball <- list()
  Cen <- vector("integer", 0)
  NotExa <- rep(TRUE, np)
  NotExa[is.na(RelSize)] <- FALSE
  Dist <- rep(1e8, np)
  BoP <- integer(np)
  nb <- 0
  
  # Create random permutation with priority to small sets
  RandPerm <- create_prioritized_permutation(RelSize)
  
  dims <- dim(Cubes)
  
  # Debug: check first few iterations
  debug_count <- 0
  nb_debug <- 0
  
  for (i in 1:np) {
    Q <- RandPerm[i]
    if (is.na(Q) || Q < 1 || Q > np || !NotExa[Q]) next
    
    # Compute set size
    # FIX: RelSize[Q] puo essere NA per punti non assegnati a segmenti
    rs_raw <- RelSize[Q]
    if (is.na(rs_raw)) rs_raw <- 128
    rs <- rs_raw / 256 * (1 - MRS) + MRS
    MaxDist <- PatchDiamMax * rs
    Radius <- MaxDist + sqrt(max(rs, 0)) * e
    # FIX: Radius deve essere positivo e finito
    if (!is.finite(Radius) || Radius <= 0) Radius <- PatchDiamMax
    # FIX (vedi conversazione): fedele all'originale N=ceil(Radius/r),
    # usando il vero r calcolato per la partizione (passato come
    # parametro) invece dell'approssimazione errata PatchDiamMax/4 -
    # verificato empiricamente che il tempo di esecuzione sui casi piu'
    # difficili (es. l'albero 47, notoriamente contaminato) e' identico
    # al comportamento preesistente, non una regressione introdotta da
    # questo fix.
    N <- max(1L, ceiling(Radius / r))
    
    # Get neighboring cubes
    cx <- CC[Q, 1]
    cy <- CC[Q, 2]
    cz <- CC[Q, 3]
    
    points <- c()
    for (dx in -N:N) {
      for (dy in -N:N) {
        for (dz in -N:N) {
          nx <- cx + dx
          ny <- cy + dy
          nz <- cz + dz
          if (nx >= 1 && nx <= dims[1] && 
              ny >= 1 && ny <= dims[2] && 
              nz >= 1 && nz <= dims[3]) {
            cell <- Partition[nx, ny, nz]
            if (!is.null(cell) && length(cell) > 0) {
              points <- c(points, cell[[1]])
            }
          }
        }
      }
    }
    
    if (length(points) == 0) next
    
    # Make sure points is valid
    if (is.list(points)) {
      points <- unlist(points)
    }
    points <- as.integer(points)
    points <- points[!is.na(points)]
    points <- points[points > 0 & points <= np]
    points <- unique(points)
    
    if (length(points) == 0) next
    
    # Compute distances
    V <- P[points, , drop = FALSE] - matrix(P[Q, ], nrow = length(points), ncol = 3, byrow = TRUE)
    dist <- rowSums(V^2)
    
    # Select points inside ball
    Inside <- dist < Radius^2
    if (sum(Inside) >= nmin) {
      ball <- points[Inside]
      d <- dist[Inside]
      core <- d < MaxDist^2
      NotExa[ball[core]] <- FALSE
      
      # Define new ball
      nb <- nb + 1
      Ball[[nb]] <- ball
      Cen[nb] <- Q
      
      # Update distances
      D <- Dist[ball]
      closer <- d < D
      ball <- ball[closer]
      Dist[ball] <- d[closer]
      BoP[ball] <- nb
    }
  }
  
  # Trim vectors
  if (nb > 0) {
    Ball <- Ball[1:nb]
    Cen <- Cen[1:nb]
  }
  
  # Build cover sets
  cover <- build_cover_structure(P, Ball, Cen, nb, BoP)
  
  return(cover)
}

#' Create prioritized permutation for variable-size cover
create_prioritized_permutation <- function(RelSize) {
  np <- length(RelSize)
  ind <- 1:np
  
  # Group by size
  small <- which(RelSize <= 32)
  medium <- which(RelSize > 32 & RelSize <= 128)
  large <- which(RelSize > 128)
  
  RandPerm <- integer(np)
  
  if (length(small) > 0) {
    RandPerm[1:length(small)] <- sample(small)
  }
  
  if (length(medium) > 0) {
    start <- length(small) + 1
    RandPerm[start:(start + length(medium) - 1)] <- sample(medium)
  }
  
  if (length(large) > 0) {
    start <- length(small) + length(medium) + 1
    RandPerm[start:np] <- sample(large)
  }
  
  return(RandPerm)
}

#' Build cover structure with neighbors
#'
#' Porting fedele della seconda meta' di cover_sets.m (sezioni "Cover
#' sets" e "Neighbors"). FIX (vedi conversazione): riceve ora il vero
#' BoP calcolato durante il loop principale di generazione (che traccia
#' quale seed e' effettivamente "piu' vicino" per ciascun punto,
#' aggiornato dinamicamente man mano che nuovi ball vengono creati) -
#' prima veniva ricostruito qui usando semplicemente Ball[[b]] (il
#' "big ball" salvato al momento della creazione, MAI aggiornato per i
#' cambi successivi di BoP), violando la proprieta' fondamentale del
#' cover ("ogni punto appartiene al massimo un cover set", come dice
#' esplicitamente il commento dell'originale) - un punto poteva
#' comparire in piu' di un cover set contemporaneamente.
#'
#' FIX (secondo problema correlato): Nei (i vicini) ora viene costruito
#' usando Ball (il "big ball" CON overlap, come nell'originale) per
#' trovare i punti candidati, poi mappandoli al loro vero cover set
#' tramite BoP - prima veniva usato erroneamente PointsInSets (i punti
#' ESCLUSIVI gia' assegnati, SENZA overlap) per la ricerca dei vicini,
#' perdendo la capacita' di trovare vicini attraverso l'area di
#' sovrapposizione tra ball adiacenti, che e' esattamente il
#' meccanismo con cui l'originale definisce la relazione di vicinanza.
#'
#' @param BoP Vettore (lunghezza np) con l'indice del cover set a cui
#'   ciascun punto e' stato EFFETTIVAMENTE assegnato (0 se nessuno)
build_cover_structure <- function(P, Ball, Cen, nb, BoP) {
  np <- nrow(P)
  
  if (nb == 0 || length(Ball) == 0) {
    return(list(
      ball = list(),
      center = integer(0),
      neighbor = list()
    ))
  }
  
  Ball <- Ball[1:min(nb, length(Ball))]
  
  # Numero di punti in ciascun ball e indice di ciascun punto nel suo
  # ball, basati sul VERO BoP (fedele: sezione "Cover sets")
  Num <- integer(nb)
  Ind <- integer(np)
  for (i in seq_len(np)) {
    if (BoP[i] > 0) {
      Num[BoP[i]] <- Num[BoP[i]] + 1L
      Ind[i] <- Num[BoP[i]]
    }
  }
  
  PointsInSets <- vector("list", nb)
  for (i in seq_len(nb)) {
    PointsInSets[[i]] <- integer(Num[i])
  }
  for (i in seq_len(np)) {
    if (BoP[i] > 0) {
      PointsInSets[[BoP[i]]][Ind[i]] <- i
    }
  }
  
  # Vicini: usa Ball (il "big ball" CON overlap), come nell'originale
  # (fedele: sezione "Neighbors")
  Nei <- vector("list", nb)
  for (i in seq_len(nb)) {
    B <- Ball[[i]]
    if (is.null(B) || length(B) == 0) {
      Nei[[i]] <- integer(0)
      next
    }
    B <- B[B > 0 & B <= np]
    if (length(B) == 0) {
      Nei[[i]] <- integer(0)
      next
    }
    I <- BoP[B] != i
    N <- B[I]
    N <- BoP[N]
    N <- N[N > 0]
    N <- unique(N)
    Nei[[i]] <- N
  }
  
  # Rende la relazione simmetrica
  for (i in seq_len(nb)) {
    N <- Nei[[i]]
    if (length(N) == 0) next
    N <- N[N > 0 & N <= nb]
    for (j in N) {
      if (!(i %in% Nei[[j]])) {
        Nei[[j]] <- c(Nei[[j]], i)
      }
    }
  }
  
  cover <- list(
    ball = PointsInSets,
    center = Cen,
    neighbor = Nei
  )
  
  return(cover)
}