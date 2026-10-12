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

#' Update tree data structure
#'
#' @param QSM QSM structure
#' @param cylinder Cylinder structure
#' @param branch Branch structure
#' @param inputs Inputs structure
#' @return Updated treedata structure
update_tree_data <- function(QSM, cylinder, branch, inputs) {
  # Define some variables from cylinder:
  treedata <- QSM$treedata
  Rad <- cylinder$radius
  Len <- cylinder$length
  Sta <- cylinder$start
  Axe <- cylinder$axis
  nc <- length(Rad)
  ind <- (1:nc)
  Trunk <- cylinder$branch == 1  # Trunk cylinders
  
  # Tree attributes from cylinders
  # Volumes, areas, lengths, branches
  treedata$TotalVolume <- 1000 * pi * sum(Rad^2 * Len)
  treedata$TrunkVolume <- 1000 * pi * sum(Rad[Trunk]^2 * Len[Trunk])
  treedata$BranchVolume <- 1000 * pi * sum(Rad[!Trunk]^2 * Len[!Trunk])
  bottom <- min(Sta[, 3])
  top_idx <- which.max(Sta[, 3])
  top <- Sta[top_idx, 3]
  if (Axe[top_idx, 3] > 0) {
    top <- top + Len[top_idx] * Axe[top_idx, 3]
  }
  treedata$TreeHeight <- top - bottom
  treedata$TrunkLength <- sum(Len[Trunk])
  treedata$BranchLength <- sum(Len[!Trunk])
  treedata$TotalLength <- treedata$TrunkLength + treedata$BranchLength
  NB <- length(branch$order) - 1  # number of branches
  treedata$NumberBranches <- NB
  BO <- max(branch$order)  # maximum branch order
  treedata$MaxBranchOrder <- BO
  treedata$TrunkArea <- 2 * pi * sum(Rad[Trunk] * Len[Trunk])
  treedata$BranchArea <- 2 * pi * sum(Rad[!Trunk] * Len[!Trunk])
  treedata$TotalArea <- 2 * pi * sum(Rad * Len)
  
  # Crown measures, Vertical profile and spreads
  result <- crown_measures(treedata, cylinder, branch)
  treedata <- result$treedata
  spreads <- result$spreads
  
  # Update triangulation information
  if (inputs$Tria) {
    treedata <- update_triangulation(QSM, treedata, cylinder)
  }
  
  # Tree Location
  treedata$location <- Sta[1, ]
  
  # Stem taper
  R <- Rad[Trunk]
  n <- length(R)
  Taper <- matrix(0, n + 1, 2)
  Taper[1, 2] <- 2 * R[1]
  Taper[2:(n + 1), 1] <- cumsum(Len[Trunk])
  Taper[2:(n + 1), 2] <- c(2 * R[2:n], 2 * R[n])
  treedata$StemTaper <- t(Taper)
  
  # Vertical profile and spreads
  treedata$VerticalProfile <- rowMeans(spreads)
  treedata$spreads <- spreads
  
  # CYLINDER DISTRIBUTIONS:
  # Wood part diameter distributions
  treedata <- cylinder_distribution(treedata, Rad, Len, Axe, "Dia")
  
  # Wood part height distributions
  treedata <- cylinder_height_distribution(treedata, Rad, Len, Sta, Axe, ind)
  
  # Wood part zenith direction distributions
  treedata <- cylinder_distribution(treedata, Rad, Len, Axe, "Zen")
  
  # Wood part azimuth direction distributions
  treedata <- cylinder_distribution(treedata, Rad, Len, Axe, "Azi")
  
  # BRANCH DISTRIBUTIONS:
  # Branch order distributions
  treedata <- branch_order_distribution(treedata, branch)
  
  # Branch diameter distributions
  treedata <- branch_distribution(treedata, branch, "Dia")
  
  # Branch height distribution
  treedata <- branch_distribution(treedata, branch, "Hei")
  
  # Branch angle distribution
  treedata <- branch_distribution(treedata, branch, "Ang")
  
  # Branch azimuth distribution
  treedata <- branch_distribution(treedata, branch, "Azi")
  
  # Branch zenith distribution
  treedata <- branch_distribution(treedata, branch, "Zen")
  
  # Change into single-format
  Names <- names(treedata)
  n <- length(Names)
  for (i in 1:n) {
    treedata[[Names[i]]] <- as.numeric(treedata[[Names[i]]])
  }
  
  if (inputs$disp == 2) {
    # Generate units for displaying the treedata
    Units <- character(n)
    m <- 23
    for (i in 1:n) {
      if (!inputs$Tria && Names[i] == "CrownVolumeAlpha") {
        m <- i
      } else if (inputs$Tria && Names[i] == "TriaTrunkLength") {
        m <- i
      }
      if (substr(Names[i], 1, 3) == "DBH") {
        Units[i] <- "m  "
      } else if (substr(Names[i], nchar(Names[i]) - 2, nchar(Names[i])) == "ume") {
        Units[i] <- "L  "
      } else if (substr(Names[i], nchar(Names[i]) - 2, nchar(Names[i])) == "ght") {
        Units[i] <- "m  "
      } else if (substr(Names[i], nchar(Names[i]) - 2, nchar(Names[i])) == "gth") {
        Units[i] <- "m  "
      } else if (substr(Names[i], 1, 3) == "vol") {
        Units[i] <- "L  "
      } else if (substr(Names[i], 1, 3) == "len") {
        Units[i] <- "m  "
      } else if (substr(Names[i], nchar(Names[i]) - 2, nchar(Names[i])) == "rea") {
        Units[i] <- "m^2"
      } else if (substr(Names[i], 1, 3) == "loc") {
        Units[i] <- "m  "
      } else if (substr(Names[i], nchar(Names[i]) - 4, nchar(Names[i])) == "aConv") {
        Units[i] <- "m^2"
      } else if (substr(Names[i], nchar(Names[i]) - 5, nchar(Names[i])) == "aAlpha") {
        Units[i] <- "m^2"
      } else if (substr(Names[i], nchar(Names[i]) - 4, nchar(Names[i])) == "eConv") {
        Units[i] <- "m^3"
      } else if (substr(Names[i], nchar(Names[i]) - 5, nchar(Names[i])) == "eAlpha") {
        Units[i] <- "m^3"
      } else if (substr(Names[i], nchar(Names[i]) - 2, nchar(Names[i])) == "Ave") {
        Units[i] <- "m  "
      } else if (substr(Names[i], nchar(Names[i]) - 2, nchar(Names[i])) == "Max") {
        Units[i] <- "m  "
      }
    }
    # Display treedata
    cat("------------\n")
    cat("  Tree attributes:\n")
    for (i in 1:m) {
      v <- change_precision(treedata[[Names[i]]])
      if (Names[i] == "DBHtri") {
        cat("  -----\n")
        cat("  Tree attributes from triangulation:\n")
      }
      cat("  ", Names[i], " = ", v, " ", Units[i], "\n")
    }
    cat("  -----\n")
  }
  
  return(treedata)
}


#' Crown measures helper function
#'
#' @param treedata Treedata structure
#' @param cylinder Cylinder structure
#' @param branch Branch structure
#' @return List with treedata and spreads
crown_measures <- function(treedata, cylinder, branch) {
  # Generate point clouds from the cylinder model
  Axe <- cylinder$axis
  Len <- cylinder$length
  Sta <- cylinder$start
  Tip <- Sta + cbind(Len * Axe[, 1], Len * Axe[, 2], Len * Axe[, 3])  # tips of the cylinders
  nc <- length(Len)
  P <- matrix(0, 5 * nc, 3)  # four mid points on the cylinder surface
  t <- 0
  for (i in 1:nc) {
    orthonormal_result <- orthonormal_vectors(Axe[i, ])
    U <- orthonormal_result$V   # FIX: orthonormal_vectors() restituisce $V/$W
    V <- orthonormal_result$W
    U <- cylinder$radius[i] * U
    if (cylinder$branch[i] == 1) {
      # For stem cylinders generate more points
      for (k in 1:4) {
        M <- Sta[i, ] + k * Len[i] / 4 * Axe[i, ]
        R_mat <- rotation_matrix(Axe[i, ], pi / 12)
        for (j in 1:12) {
          if (j > 1) {
            U <- R_mat %*% U
          }
          t <- t + 1
          P[t, ] <- M + as.numeric(U)
        }
      }
    } else {
      M <- Sta[i, ] + Len[i] / 2 * Axe[i, ]
      R_mat <- rotation_matrix(Axe[i, ], pi / 4)
      for (j in 1:4) {
        if (j > 1) {
          U <- R_mat %*% U
        }
        t <- t + 1
        P[t, ] <- M + as.numeric(U)
      }
    }
  }
  P <- P[1:t, ]
  P <- rbind(P, Sta, Tip)
  P <- unique(P, MARGIN = 1)
  
  # Vertical profiles (layer diameters/spreads), mean:
  bot <- min(P[, 3])
  top <- max(P[, 3])
  Hei <- top - bot
  if (Hei > 10) {
    m <- 20
  } else if (Hei > 2) {
    m <- 10
  } else {
    m <- 5
  }
  spreads <- matrix(0, m, 18)
  for (j in 1:m) {
    I <- P[, 3] >= bot + (j - 1) * Hei / m & P[, 3] < bot + j * Hei / m
    X <- unique(P[I, , drop = FALSE], MARGIN = 1)
    if (nrow(X) > 5) {
      K <- chull(X[, 1], X[, 2])
      # compute center of gravity for the convex hull
      n <- length(K)
      x <- X[K, 1]
      y <- X[K, 2]
      CX <- sum((x[1:(n - 1)] + x[2:n]) * (x[1:(n - 1)] * y[2:n] - x[2:n] * y[1:(n - 1)])) / 6 / sum(chull(X[, 1], X[, 2]))
      CY <- sum((y[1:(n - 1)] + y[2:n]) * (x[1:(n - 1)] * y[2:n] - x[2:n] * y[1:(n - 1)])) / 6 / sum(chull(X[, 1], X[, 2]))
      
      V <- mat_vec_subtraction(X[, 1:2, drop = FALSE], c(CX, CY))
      ang <- atan2(V[, 2], V[, 1]) + pi
      ang_sorted <- sort(ang, index.return = TRUE)
      ang <- ang_sorted$x
      I <- ang_sorted$ix
      L <- sqrt(rowSums(V^2))
      L <- L[I]
      for (i in 1:18) {
        I1 <- ang >= (i - 1) * pi / 18 & ang < i * pi / 18
        if (any(I1)) {
          L1 <- max(L[I1])
        } else {
          L1 <- 0
        }
        J1 <- ang >= (i - 1) * pi / 18 + pi & ang < i * pi / 18 + pi
        if (any(J1)) {
          L2 <- max(L[J1])
        } else {
          L2 <- 0
        }
        spreads[j, i] <- L1 + L2
      }
    }
  }
  
  # Crown diameters (spreads), mean and maximum:
  X <- unique(P[, 1:2, drop = FALSE], MARGIN = 1)
  K <- chull(X[, 1], X[, 2])
  A <- sum(chull(X[, 1], X[, 2]))
  n <- length(K)
  x <- X[K, 1]
  y <- X[K, 2]
  CX <- sum((x[1:(n - 1)] + x[2:n]) * (x[1:(n - 1)] * y[2:n] - x[2:n] * y[1:(n - 1)])) / 6 / A
  CY <- sum((y[1:(n - 1)] + y[2:n]) * (x[1:(n - 1)] * y[2:n] - x[2:n] * y[1:(n - 1)])) / 6 / A
  V <- mat_vec_subtraction(Tip[, 1:2, drop = FALSE], c(CX, CY))
  ang <- atan2(V[, 2], V[, 1]) + pi
  ang_sorted <- sort(ang, index.return = TRUE)
  ang <- ang_sorted$x
  I <- ang_sorted$ix
  L <- sqrt(rowSums(V^2))
  L <- L[I]
  S <- numeric(18)
  for (i in 1:18) {
    I1 <- ang >= (i - 1) * pi / 18 & ang < i * pi / 18
    if (any(I1)) {
      L1 <- max(L[I1])
    } else {
      L1 <- 0
    }
    J1 <- ang >= (i - 1) * pi / 18 + pi & ang < i * pi / 18 + pi
    if (any(J1)) {
      L2 <- max(L[J1])
    } else {
      L2 <- 0
    }
    S[i] <- L1 + L2
  }
  treedata$CrownDiamAve <- mean(S)
  MaxDiam <- 0
  for (i in 1:n) {
    V <- mat_vec_subtraction(cbind(x, y), c(x[i], y[i]))
    L <- max(sqrt(rowSums(V^2)))
    if (L > MaxDiam) {
      MaxDiam <- L
    }
  }
  treedata$CrownDiamMax <- MaxDiam
  
  # Crown areas from convex hull and alpha shape:
  treedata$CrownAreaConv <- A
  alp <- max(0.5, treedata$CrownDiamAve / 10)
  # Alpha shape approximation
  treedata$CrownAreaAlpha <- A  # Simplified: using convex hull area
  
  # Crown base
  # Define first major branch
  dbh <- treedata$DBHcyl
  nb <- length(branch$order)
  HL <- numeric(nb)  # horizontal reach
  branches1 <- which(branch$order == 1)
  nb <- length(branches1)
  nc <- nrow(Sta)
  ind <- (1:nc)
  for (i in 1:nb) {
    C <- ind[cylinder$branch == branches1[i]]
    if (length(C) > 0) {
      base <- Sta[C[1], ]
      C <- C[length(C)]
      tip <- Sta[C, ] + Len[C] * Axe[C, ]
      V <- tip[1:2] - base[1:2]
      HL[branches1[i]] <- sqrt(V %*% V) / dbh * 2
    }
  }
  M <- min(10, median(HL[HL > 0]))
  
  # Sort the branches according to their heights
  Hei <- branch$height[branches1]
  Hei_sorted <- sort(Hei, index.return = TRUE)
  Hei <- Hei_sorted$x
  branches1 <- branches1[Hei_sorted$ix]
  
  # Search the first/lowest branch
  d <- min(0.05, 0.05 * dbh)
  b <- 0
  if (nb > 1) {
    i <- 1
    while (i < nb) {
      i <- i + 1
      if (branch$diameter[branches1[i]] > d && HL[branches1[i]] > M) {
        b <- branches1[i]
        i <- nb + 2
      }
    }
    if (i == nb + 1 && nb > 1) {
      b <- branches1[1]
    }
  }
  
  if (b > 0) {
    # Search all the children of the first major branch
    nb <- length(branch$parent)
    Ind <- (1:nb)
    chi <- Ind[branch$parent == b]
    B <- b
    while (length(chi) > 0) {
      B <- c(B, chi)
      n <- length(chi)
      C <- vector("list", n)
      for (i in 1:n) {
        C[[i]] <- Ind[branch$parent == chi[i]]
      }
      chi <- unlist(C)
    }
    
    # Define crown base height from the ground
    BaseHeight <- max(Sta[, 3])  # Height of the crown base
    for (i in 1:length(B)) {
      C <- ind[cylinder$branch == B[i]]
      ht <- min(Tip[C, 3])
      hb <- min(Sta[C, 3])
      h <- min(hb, ht)
      if (h < BaseHeight) {
        BaseHeight <- h
      }
    }
    treedata$CrownBaseHeight <- BaseHeight - Sta[1, 3]
    
    # Crown length and ratio
    treedata$CrownLength <- treedata$TreeHeight - treedata$CrownBaseHeight
    treedata$CrownRatio <- treedata$CrownLength / treedata$TreeHeight
    
    # Crown volume from convex hull
    I <- P[, 3] >= BaseHeight
    X <- P[I, , drop = FALSE]
    # Simplified volume calculation
    treedata$CrownVolumeConv <- 0
    treedata$CrownVolumeAlpha <- 0
    
  } else {
    # No branches
    treedata$CrownBaseHeight <- treedata$TreeHeight
    treedata$CrownLength <- 0
    treedata$CrownRatio <- 0
    treedata$CrownVolumeConv <- 0
    treedata$CrownVolumeAlpha <- 0
  }
  
  return(list(treedata = treedata, spreads = spreads))
}


#' Update triangulation helper
#'
#' @param QSM QSM structure
#' @param treedata Treedata structure
#' @param cylinder Cylinder structure
#' @return Updated treedata
update_triangulation <- function(QSM, treedata, cylinder) {
  # Update the mixed results
  if (!is.null(QSM$triangulation)) {
    CylInd <- QSM$triangulation$cylind
    Rad <- cylinder$radius
    Len <- cylinder$length
    # Determine the new stem cylinder
    nc <- length(Rad)
    ind <- (1:nc)
    ind <- ind[cylinder$branch == 1]  # cylinders in the stem
    S <- QSM$cylinder$start[CylInd, ]  # The place where the triangulation stops
    V <- cylinder$start[ind, , drop = FALSE] - S
    d <- sqrt(rowSums(V^2))
    d_min <- min(d)
    I <- which.min(d)
    V <- V[I, ]
    CylInd <- ind[I]  # The new cylinder closest to the correct place
    
    if (d_min < 0.01) {
      TrunkVolMix <- treedata$TrunkVolume - 
        1000 * pi * sum(Rad[1:(CylInd - 1)]^2 * Len[1:(CylInd - 1)]) + QSM$triangulation$volume
      TrunkAreaMix <- treedata$TrunkArea - 
        2 * pi * sum(Rad[1:(CylInd - 1)] * Len[1:(CylInd - 1)]) + QSM$triangulation$SideArea
    } else {
      # Select the following cylinder
      h <- V %*% cylinder$axis[CylInd, ]
      if (h < 0) {
        CylInd <- CylInd + 1
        V <- cylinder$start[CylInd, ] - S
        h <- V %*% cylinder$axis[CylInd, ]
      }
      Len[CylInd - 1] <- Len[CylInd - 1] - h
      
      TrunkVolMix <- treedata$TrunkVolume - 
        1000 * pi * sum(Rad[1:(CylInd - 1)]^2 * Len[1:(CylInd - 1)]) + QSM$triangulation$volume
      TrunkAreaMix <- treedata$TrunkArea - 
        2 * pi * sum(Rad[1:(CylInd - 1)] * Len[1:(CylInd - 1)]) + QSM$triangulation$SideArea
    }
    treedata$MixTrunkVolume <- TrunkVolMix
    treedata$MixTotalVolume <- TrunkVolMix + treedata$BranchVolume
    treedata$MixTrunkArea <- TrunkAreaMix
    treedata$MixTotalArea <- TrunkAreaMix + treedata$BranchArea
  }
  return(treedata)
}


#' Cylinder distribution helper
#'
#' @param treedata Treedata structure
#' @param Rad Radius vector
#' @param Len Length vector
#' @param Axe Axis matrix
#' @param dist Distribution type ("Dia", "Zen", "Azi")
#' @return Updated treedata
cylinder_distribution <- function(treedata, Rad, Len, Axe, dist) {
  if (dist == "Dia") {
    Par <- Rad
    n <- ceiling(max(200 * Rad))
    a <- 0.005  # diameter in 1 cm classes
  } else if (dist == "Zen") {
    Par <- 180 / pi * acos(Axe[, 3])
    n <- 18
    a <- 10  # zenith direction in 10 degree angle classes
  } else if (dist == "Azi") {
    Par <- 180 / pi * atan2(Axe[, 2], Axe[, 1]) + 180
    n <- 36
    a <- 10  # azimuth direction in 10 degree angle classes
  }
  
  CylDist <- matrix(0, 3, n)
  for (i in 1:n) {
    K <- Par >= (i - 1) * a & Par < i * a
    CylDist[1, i] <- 1000 * pi * sum(Rad[K]^2 * Len[K])  # volumes in litres
    CylDist[2, i] <- 2 * pi * sum(Rad[K] * Len[K])  # areas in litres
    CylDist[3, i] <- sum(Len[K])  # lengths in meters
  }
  treedata[[paste0("VolCyl", dist)]] <- CylDist[1, ]
  treedata[[paste0("AreCyl", dist)]] <- CylDist[2, ]
  treedata[[paste0("LenCyl", dist)]] <- CylDist[3, ]
  return(treedata)
}


#' Cylinder height distribution helper
#'
#' @param treedata Treedata structure
#' @param Rad Radius vector
#' @param Len Length vector
#' @param Sta Start matrix
#' @param Axe Axis matrix
#' @param ind Index vector
#' @return Updated treedata
cylinder_height_distribution <- function(treedata, Rad, Len, Sta, Axe, ind) {
  MaxHei <- ceiling(treedata$TreeHeight)
  treedata$VolCylHei <- numeric(MaxHei)
  treedata$AreCylHei <- numeric(MaxHei)
  treedata$LenCylHei <- numeric(MaxHei)
  End <- Sta + cbind(Len * Axe[, 1], Len * Axe[, 2], Len * Axe[, 3])
  bot <- min(Sta[, 3])
  B <- Sta[, 3] - bot
  T <- End[, 3] - bot
  for (j in 1:MaxHei) {
    I1 <- B >= (j - 2) & B < (j - 1)  # base below this bin
    J1 <- B >= (j - 1) & B < j  # base in this bin
    K1 <- B >= j & B < (j + 1)  # base above this bin
    I2 <- T >= (j - 2) & T < (j - 1)  # top below this bin
    J2 <- T >= (j - 1) & T < j  # top in this bin
    K2 <- T >= j & T < (j + 1)  # top above this bin
    C1 <- ind[J1 & J2]  # base and top in this bin
    C2 <- ind[J1 & K2]  # base in this bin, top above
    C3 <- ind[J1 & I2]  # base in this bin, top below
    C4 <- ind[I1 & J2]  # base in bin below, top in this
    C5 <- ind[K1 & J2]  # base in bin above, top in this
    
    v1 <- 1000 * pi * sum(Rad[C1]^2 * Len[C1])
    a1 <- 2 * pi * sum(Rad[C1] * Len[C1])
    l1 <- sum(Len[C1])
    
    if (length(C2) > 0) {
      r2 <- (j - B[C2]) / (T[C2] - B[C2])  # relative portion in this bin
      v2 <- 1000 * pi * sum(Rad[C2]^2 * Len[C2] * r2)
      a2 <- 2 * pi * sum(Rad[C2] * Len[C2] * r2)
      l2 <- sum(Len[C2] * r2)
    } else {
      v2 <- 0
      a2 <- 0
      l2 <- 0
    }
    
    if (length(C3) > 0) {
      r3 <- (B[C3] - j + 1) / (B[C3] - T[C3])  # relative portion in this bin
      v3 <- 1000 * pi * sum(Rad[C3]^2 * Len[C3] * r3)
      a3 <- 2 * pi * sum(Rad[C3] * Len[C3] * r3)
      l3 <- sum(Len[C3] * r3)
    } else {
      v3 <- 0
      a3 <- 0
      l3 <- 0
    }
    
    if (length(C4) > 0) {
      r4 <- (T[C4] - j + 1) / (T[C4] - B[C4])  # relative portion in this bin
      v4 <- 1000 * pi * sum(Rad[C4]^2 * Len[C4] * r4)
      a4 <- 2 * pi * sum(Rad[C4] * Len[C4] * r4)
      l4 <- sum(Len[C4] * r4)
    } else {
      v4 <- 0
      a4 <- 0
      l4 <- 0
    }
    
    if (length(C5) > 0) {
      r5 <- (j - T[C5]) / (B[C5] - T[C5])  # relative portion in this bin
      v5 <- 1000 * pi * sum(Rad[C5]^2 * Len[C5] * r5)
      a5 <- 2 * pi * sum(Rad[C5] * Len[C5] * r5)
      l5 <- sum(Len[C5] * r5)
    } else {
      v5 <- 0
      a5 <- 0
      l5 <- 0
    }
    
    treedata$VolCylHei[j] <- v1 + v2 + v3 + v4 + v5
    treedata$AreCylHei[j] <- a1 + a2 + a3 + a4 + a5
    treedata$LenCylHei[j] <- l1 + l2 + l3 + l4 + l5
  }
  return(treedata)
}


#' Branch distribution helper
#'
#' @param treedata Treedata structure
#' @param branch Branch structure
#' @param dist Distribution type ("Dia", "Hei", "Ang", "Zen", "Azi")
#' @return Updated treedata
branch_distribution <- function(treedata, branch, dist) {
  BOrd <- branch$order[-1]
  BVol <- branch$volume[-1]
  BAre <- branch$area[-1]
  BLen <- branch$length[-1]
  
  if (dist == "Dia") {
    Par <- branch$diameter[-1]
    n <- ceiling(max(100 * Par))
    a <- 0.005  # diameter in 1 cm classes
  } else if (dist == "Hei") {
    Par <- branch$height[-1]
    n <- ceiling(treedata$TreeHeight)
    a <- 1  # height in 1 m classes
  } else if (dist == "Ang") {
    Par <- branch$angle[-1]
    n <- 18
    a <- 10  # angle in 10 degree classes
  } else if (dist == "Zen") {
    Par <- branch$zenith[-1]
    n <- 18
    a <- 10  # zenith direction in 10 degree angle classes
  } else if (dist == "Azi") {
    Par <- branch$azimuth[-1] + 180
    n <- 36
    a <- 10  # azimuth direction in 10 degree angle classes
  }
  
  BranchDist <- matrix(0, 8, n)
  for (i in 1:n) {
    I <- Par >= (i - 1) * a & Par < i * a
    BranchDist[1, i] <- sum(BVol[I])  # volume (all branches)
    BranchDist[2, i] <- sum(BVol[I & BOrd == 1])  # volume (1st-branches)
    BranchDist[3, i] <- sum(BAre[I])  # area (all branches)
    BranchDist[4, i] <- sum(BAre[I & BOrd == 1])  # area (1st-branches)
    BranchDist[5, i] <- sum(BLen[I])  # length (all branches)
    BranchDist[6, i] <- sum(BLen[I & BOrd == 1])  # length (1st-branches)
    BranchDist[7, i] <- sum(I)  # number (all branches)
    BranchDist[8, i] <- sum(I & BOrd == 1)  # number (1st-branches)
  }
  treedata[[paste0("VolBranch", dist)]] <- BranchDist[1, ]
  treedata[[paste0("VolBranch1", dist)]] <- BranchDist[2, ]
  treedata[[paste0("AreBranch", dist)]] <- BranchDist[3, ]
  treedata[[paste0("AreBranch1", dist)]] <- BranchDist[4, ]
  treedata[[paste0("LenBranch", dist)]] <- BranchDist[5, ]
  treedata[[paste0("LenBranch1", dist)]] <- BranchDist[6, ]
  treedata[[paste0("NumBranch", dist)]] <- BranchDist[7, ]
  treedata[[paste0("NumBranch1", dist)]] <- BranchDist[8, ]
  return(treedata)
}


#' Branch order distribution helper
#'
#' @param treedata Treedata structure
#' @param branch Branch structure
#' @return Updated treedata
branch_order_distribution <- function(treedata, branch) {
  BO <- max(branch$order)
  BranchOrdDist <- matrix(0, max(1, BO), 4)
  for (i in 1:max(1, BO)) {
    I <- branch$order == i
    BranchOrdDist[i, 1] <- sum(branch$volume[I])  # volumes
    BranchOrdDist[i, 2] <- sum(branch$area[I])  # areas
    BranchOrdDist[i, 3] <- sum(branch$length[I])  # lengths
    BranchOrdDist[i, 4] <- sum(I)  # number of ith-order branches
  }
  treedata$VolBranchOrd <- BranchOrdDist[, 1]
  treedata$AreBranchOrd <- BranchOrdDist[, 2]
  treedata$LenBranchOrd <- BranchOrdDist[, 3]
  treedata$NumBranchOrd <- BranchOrdDist[, 4]
  return(treedata)
}