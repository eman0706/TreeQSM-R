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

#' Branches
#'
#' Porting fedele di branches.m (Version 3.0.0, Copyright (C) 2013-2022
#' Pasi Raumonen) - sostituisce la versione precedente (mai verificata
#' contro il vero sorgente in questa sessione), che aveva due
#' divergenze reali:
#' 1) Mancava completamente il campo "zenith" (angolo dalla verticale)
#'    nell'output.
#' 2) La logica di parentela tra branch era ALGORITMICAMENTE DIVERSA:
#'    la versione precedente guardava il GENITORE di ciascun cilindro
#'    (cylinder$parent[i]) per dedurre il branch genitore - il vero
#'    algoritmo guarda invece i FIGLI di ciascun cilindro nel branch,
#'    escludendo esplicitamente cylinder.extension (la naturale
#'    continuazione dello stesso ramo fisico, non un vero nuovo ramo
#'    che si dirama) - due algoritmi diversi che possono divergere,
#'    specialmente quando l'estensione conta.
#'
#' @param cylinder Output di cylinders() (radius, length, axis, start,
#'   branch, BranchOrder, parent, extension, added)
#' @return Lista "branch" con order, parent, volume, area, length,
#'   angle, height, azimuth, diameter, zenith, children, n_branch
#' @export
branches <- function(cylinder) {
  Rad <- cylinder$radius
  Len <- cylinder$length
  Axe <- cylinder$axis
  
  nc <- length(Rad)
  if (nc == 0 || length(cylinder$branch) == 0) {
    return(list(
      order = integer(0), parent = integer(0), volume = numeric(0),
      area = numeric(0), length = numeric(0), angle = numeric(0),
      height = numeric(0), azimuth = numeric(0), diameter = numeric(0),
      zenith = numeric(0), children = list(), n_branch = 0L
    ))
  }
  ns <- max(cylinder$branch)
  if (ns == 0) {
    return(list(
      order = integer(0), parent = integer(0), volume = numeric(0),
      area = numeric(0), length = numeric(0), angle = numeric(0),
      height = numeric(0), azimuth = numeric(0), diameter = numeric(0),
      zenith = numeric(0), children = list(), n_branch = 0L
    ))
  }
  
  ind <- seq_len(nc)
  BData <- matrix(0, ns, 9)
  CiB <- vector("list", ns)
  
  for (i in seq_len(ns)) {
    C <- ind[cylinder$branch == i]
    CiB[[i]] <- C
    
    if (length(C) > 0) {
      BData[i, 1] <- cylinder$BranchOrder[C[1]]
      BData[i, 2] <- 2 * Rad[C[1]]
      BData[i, 3] <- 1000 * pi * sum(Len[C] * Rad[C]^2)
      BData[i, 4] <- 2 * pi * sum(Len[C] * Rad[C])
      BData[i, 5] <- sum(Len[C])
      
      if (isTRUE(cylinder$added[C[1]]) && length(C) > 1) {
        FC <- C[2]
        PC <- cylinder$parent[C[1]]
      } else {
        FC <- C[1]
        PC <- cylinder$parent[FC]
      }
      if (!is.na(PC) && PC > 0) {
        cos_angle <- sum(Axe[FC, ] * Axe[PC, ])
        BData[i, 6] <- 180 / pi * acos(max(-1, min(1, cos_angle)))
      }
      
      BData[i, 7] <- cylinder$start[C[1], 3] - cylinder$start[1, 3]
      BData[i, 8] <- 180 / pi * atan2(Axe[C[1], 2], Axe[C[1], 1])
      BData[i, 9] <- 180 / pi * acos(max(-1, min(1, Axe[C[1], 3])))
    }
  }
  
  ## Struttura di ramificazione (topologia, relazione genitore-figlio)
  branch_order <- as.integer(BData[, 1])
  BPar <- integer(ns)
  Chi <- vector("list", nc)
  
  for (i in seq_len(nc)) {
    c_idx <- ind[!is.na(cylinder$parent) & cylinder$parent == i]
    ext_i <- cylinder$extension[i]
    if (!is.na(ext_i) && ext_i != 0) {
      c_idx <- c_idx[c_idx != ext_i]
    }
    Chi[[i]] <- c_idx
  }
  for (i in seq_len(ns)) {
    C <- CiB[[i]]
    if (length(C) == 0) next
    ChildCyls <- unique(unlist(Chi[C], use.names = FALSE))
    if (length(ChildCyls) == 0) next
    CB <- unique(cylinder$branch[ChildCyls])
    CB <- CB[!is.na(CB) & CB >= 1 & CB <= ns]
    if (length(CB) > 0) BPar[CB] <- i
  }
  
  list(
    order    = branch_order,
    parent   = as.integer(BPar),
    volume   = as.numeric(BData[, 3]),
    area     = as.numeric(BData[, 4]),
    length   = as.numeric(BData[, 5]),
    angle    = as.numeric(BData[, 6]),
    height   = as.numeric(BData[, 7]),
    azimuth  = as.numeric(BData[, 8]),
    diameter = as.numeric(BData[, 2]),
    zenith   = as.numeric(BData[, 9]),
    children = Chi,
    n_branch = ns
  )
}