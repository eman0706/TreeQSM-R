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

#' Segments
#'
#' Porting fedele di segments.m (Version 2.10, Copyright (C) 2013-2017
#' Pasi Raumonen). Segmenta l'albero in rami e registra le relazioni
#' padre-figlio. Le biforcazioni sono riconosciute studiando la
#' connettivita' di una regione "study" che avanza lungo l'albero. In
#' caso di componenti multiple connesse in "study", le componenti sono
#' classificate come continuazione o come rami.
#'
#' @param cover Lista con neighbor (e altri campi di cover_sets)
#' @param Base Base dell'albero (indici cover set)
#' @param Forb Cover set non parte dell'albero
#' @return Lista con segments, ParentSegment, ChildSegment
#' @export
segments <- function(cover, Base, Forb) {
  Nei <- cover$neighbor
  nb  <- length(Nei)
  
  a    <- max(200000L, as.integer(nb / 100L))
  SBas <- vector("list", a)
  Segs <- vector("list", a)
  SPar <- matrix(0L, a, 2)
  SChi <- vector("list", a)
  NChi <- integer(a)
  
  Fal <- rep(FALSE, nb)
  
  s <- 1L
  b <- s
  
  SBas[[s]] <- Base
  Seg      <- vector("list", 2000L)
  Seg[[1]] <- Base
  
  ForbAll       <- Fal
  ForbAll[Forb] <- TRUE
  ForbAll[Base] <- TRUE
  Forb_cur      <- ForbAll
  
  Continue <- TRUE
  NewSeg   <- TRUE
  nl       <- 1L
  ns       <- 6L
  
  while (Continue && (b < nb)) {
    Forb_cur[Seg[[nl]]] <- TRUE
    
    Cut     <- define_cut(Nei, Seg[[nl]], Forb_cur, Fal)
    CutSize <- length(Cut)
    
    if (NewSeg) {
      NewSeg <- FALSE
      ns     <- min(CutSize, 6L)
    }
    
    nc <- 0L
    if (CutSize > 0L) {
      cc       <- cut_components(Nei, Cut, CutSize, Fal, Fal)
      CutComps <- cc$Components
      nc       <- length(CutComps)
      if (nc > 1L) {
        sc         <- study_components(Nei, ns, Cut, CutComps, Forb_cur, Fal, Fal)
        StudyComps <- sc$Components
        Bases_sc   <- sc$Bases
        CompSize   <- sc$CompSize
        Cont       <- sc$Cont
        BaseSize   <- sc$BaseSize
        nc         <- length(Cont)
      }
    }
    
    if (nc == 1L) {
      nl        <- nl + 1L
      Seg[[nl]] <- Cut
      
    } else if (nc > 1L) {
      Class <- component_classification(CompSize, Cont, BaseSize, CutSize)
      
      for (i in seq_len(nc)) {
        if (Class[i] == 1L) {
          Base_i <- Bases_sc[[i]]
          ForbAll[Base_i] <- TRUE
          Forb_cur[StudyComps[[i]]] <- TRUE
          J   <- Forb_cur[Cut]
          Cut <- Cut[!J]
          b <- b + 1L
          SBas[[b]] <- Base_i
          SPar[b, ] <- c(s, nl)
          NChi[s]   <- NChi[s] + 1L
          SChi[[s]] <- c(SChi[[s]], b)
        }
      }
      
      if (length(Cut) == 0L) {
        Segs[[s]] <- Seg[seq_len(nl)]
        S <- unlist(Seg[seq_len(nl)])
        ForbAll[S] <- TRUE
        
        if (s < b) {
          s        <- s + 1L
          Seg[[1]] <- SBas[[s]]
          Forb_cur <- ForbAll
          NewSeg   <- TRUE
          nl       <- 1L
        } else {
          Continue <- FALSE
        }
      } else {
        nl        <- nl + 1L
        Seg[[nl]] <- Cut
      }
      
    } else {
      Segs[[s]] <- Seg[seq_len(nl)]
      S <- unlist(Seg[seq_len(nl)])
      ForbAll[S] <- TRUE
      
      if (s < b) {
        s        <- s + 1L
        Seg[[1]] <- SBas[[s]]
        Forb_cur <- ForbAll
        NewSeg   <- TRUE
        nl       <- 1L
      } else {
        Continue <- FALSE
      }
    }
  }
  
  Segs <- Segs[seq_len(b)]
  SPar <- SPar[seq_len(b), , drop = FALSE]
  
  SChi_out <- vector("list", b)
  for (i in seq_len(b)) {
    if (NChi[i] > 0L && !is.null(SChi[[i]])) {
      SChi_out[[i]] <- as.integer(SChi[[i]])
    } else {
      SChi_out[[i]] <- integer(0)
    }
    if (!is.null(Segs[[i]])) {
      for (j in seq_along(Segs[[i]])) {
        Segs[[i]][[j]] <- as.integer(Segs[[i]][[j]])
      }
    }
  }
  
  n_layers <- sapply(Segs, function(x) if (!is.null(x)) length(x) else 0L)
  cat(sprintf("  [DIAG] segment: %d segmenti (porting fedele segments.m)\n", b))
  cat(sprintf("  [DIAG] layer/seg: min=%d median=%d max=%d\n",
              as.integer(min(n_layers)), as.integer(median(n_layers)), as.integer(max(n_layers))))
  
  list(
    segments      = Segs,
    ParentSegment = SPar,
    ChildSegment  = SChi_out
  )
}

#' @export
define_cut <- function(Nei, CutPre, Forb, Fal) {
  Cut <- unlist(Nei[CutPre])
  if (length(Cut) == 0L) return(integer(0))
  Cut <- unique(Cut)
  I <- Forb[Cut]
  Cut[!I]
}

#' Componenti connesse del taglio (Cut)
#' @export
cut_components <- function(Nei, Cut, CutSize, Fal, False) {
  if (CutSize == 1L) {
    return(list(Components = list(Cut), CompSize = 1L))
    
  } else if (CutSize == 2L) {
    I <- Nei[[Cut[1]]] == Cut[2]
    if (any(I)) {
      return(list(Components = list(Cut), CompSize = 1L))
    } else {
      return(list(Components = list(Cut[1], Cut[2]), CompSize = c(1L, 1L)))
    }
    
  } else if (CutSize == 3L) {
    I <- Nei[[Cut[1]]] == Cut[2]
    J <- Nei[[Cut[1]]] == Cut[3]
    K <- Nei[[Cut[2]]] == Cut[3]
    nI <- any(I); nJ <- any(J); nK <- any(K)
    if ((nI + nJ + nK) >= 2) {
      return(list(Components = list(Cut), CompSize = 1L))
    } else if (nI) {
      return(list(Components = list(Cut[1:2], Cut[3]), CompSize = c(2L, 1L)))
    } else if (nJ) {
      return(list(Components = list(Cut[c(1, 3)], Cut[2]), CompSize = c(2L, 1L)))
    } else if (nK) {
      return(list(Components = list(Cut[2:3], Cut[1]), CompSize = c(2L, 1L)))
    } else {
      return(list(Components = list(Cut[1], Cut[2], Cut[3]), CompSize = c(1L, 1L, 1L)))
    }
    
  } else {
    Components <- vector("list", CutSize)
    CompSize   <- integer(CutSize)
    Comp       <- integer(CutSize)
    Fal[Cut]   <- TRUE
    nc <- 0L
    m  <- Cut[1]
    i  <- 0L
    while (i < CutSize) {
      Added <- Nei[[m]]
      I     <- Fal[Added]
      Added <- Added[I]
      a     <- length(Added)
      Comp[1] <- m
      Fal[m]  <- FALSE
      t <- 1L
      while (a > 0L) {
        Comp[(t + 1):(t + a)] <- Added
        Fal[Added] <- FALSE
        t <- t + a
        Ext <- unique(unlist(Nei[Added]))
        I <- Fal[Ext]
        Added <- Ext[I]
        a <- length(Added)
      }
      i  <- i + t
      nc <- nc + 1L
      Components[[nc]] <- Comp[seq_len(t)]
      CompSize[nc]      <- t
      if (i < CutSize) {
        J <- Fal[Cut]
        m <- Cut[J][1]
      }
    }
    Components <- Components[seq_len(nc)]
    CompSize   <- CompSize[seq_len(nc)]
    list(Components = Components, CompSize = CompSize)
  }
}

#' Componenti della regione "study" (Cut + layer successivi lungo l'albero)
#' @export
study_components <- function(Nei, ns, Cut, CutComps, Forb, Fal, False) {
  Study     <- vector("list", ns)
  StudySize <- integer(ns)
  Study[[1]]    <- Cut
  StudySize[1]  <- length(Cut)
  
  if (ns >= 2L) {
    N <- Cut
    i <- 1L
    while (i < ns) {
      Forb[N] <- TRUE
      N <- unique(unlist(Nei[N]))
      I <- Forb[N]
      N <- N[!I]
      if (length(N) > 0L) {
        i <- i + 1L
        Study[[i]]   <- N
        StudySize[i] <- length(N)
      } else {
        Study     <- Study[seq_len(i)]
        StudySize <- StudySize[seq_len(i)]
        i <- ns + 1L
      }
    }
  }
  
  ns2       <- length(StudySize)
  studysize <- sum(StudySize)
  study     <- unlist(Study)
  
  nc <- length(CutComps)
  i <- 1L; j <- 0L; k <- 0L
  Fal[study] <- TRUE
  Components <- vector("list", nc)
  CompSize   <- integer(nc)
  Comp       <- integer(studysize)
  
  while (i <= nc) {
    C <- CutComps[[i]]
    while (j < studysize) {
      a <- length(C)
      Comp[seq_len(a)] <- C
      Fal[C] <- FALSE
      if (a > 1L) {
        Add <- unique(unlist(Nei[C]))
      } else {
        Add <- Nei[[C]]
      }
      t <- a
      I <- Fal[Add]
      Add <- Add[I]
      a <- length(Add)
      while (a > 0L) {
        Comp[(t + 1):(t + a)] <- Add
        Fal[Add] <- FALSE
        t <- t + a
        Add <- unique(unlist(Nei[Add]))
        I <- Fal[Add]
        Add <- Add[I]
        a <- length(Add)
      }
      j <- j + t
      k <- k + 1L
      Components[[k]] <- Comp[seq_len(t)]
      CompSize[k]      <- t
      
      if (j < studysize) {
        C <- integer(0)
        while (i < nc && length(C) == 0L) {
          i <- i + 1L
          C <- CutComps[[i]]
          J <- Fal[C]
          C <- C[J]
        }
        if (i == nc && length(C) == 0L) {
          j <- studysize
          i <- nc + 1L
        }
      } else {
        i <- nc + 1L
      }
    }
  }
  Components <- Components[seq_len(k)]
  CompSize   <- CompSize[seq_len(k)]
  
  Cont     <- rep(TRUE, k)
  BaseSize <- integer(k)
  Bases    <- vector("list", k)
  
  if (k > 1L) {
    Forb[study] <- TRUE
    Fal[study]  <- FALSE
    Fal[Study[[1]]] <- TRUE
    for (i in seq_len(k)) {
      Set <- unique(c(Components[[i]], Study[[1]]))
      False[Components[[i]]] <- TRUE
      I <- False[Set] & Fal[Set]
      False[Components[[i]]] <- FALSE
      Set <- Set[I]
      Bases[[i]]   <- Set
      BaseSize[i]  <- length(Set)
    }
    Fal[Study[[1]]]   <- FALSE
    Fal[Study[[ns2]]] <- TRUE
    Forb[study] <- TRUE
    for (i in seq_len(k)) {
      Set <- unique(c(Components[[i]], Study[[ns2]]))
      False[Components[[i]]] <- TRUE
      I <- False[Set] & Fal[Set]
      False[Components[[i]]] <- FALSE
      Set <- Set[I]
      if (length(Set) > 0L) {
        N <- unique(unlist(Nei[Set]))
        I <- Forb[N]
        N <- N[!I]
        if (length(N) == 0L) Cont[i] <- FALSE
      } else {
        Cont[i] <- FALSE
      }
    }
  }
  
  list(Components = Components, Bases = Bases, CompSize = CompSize,
       Cont = Cont, BaseSize = BaseSize)
}

#' Classifica le componenti della regione "study": 0 = continuazione, 1 = ramo
#' @export
component_classification <- function(CompSize, Cont, BaseSize, CutSize) {
  nc        <- length(CompSize)
  StudySize <- sum(CompSize)
  Class     <- rep(1L, nc)
  ContiComp <- 0L
  
  for (i in seq_len(nc)) {
    if (BaseSize[i] == CompSize[i] && !Cont[i]) {
      Class[i] <- 0L
    } else if (BaseSize[i] == 1L && CompSize[i] <= 2L && !Cont[i]) {
      Class[i] <- 0L
    } else if (BaseSize[i] / CutSize < 0.05 && 2L * BaseSize[i] >= CompSize[i] && !Cont[i]) {
      Class[i] <- 0L
    } else if (CompSize[i] <= 3L && !Cont[i]) {
      Class[i] <- 0L
    } else if (BaseSize[i] / CutSize >= 0.7 || CompSize[i] >= 0.7 * StudySize) {
      Class[i]  <- 0L
      ContiComp <- i
    }
  }
  
  branches <- which(Class == 1L)
  if (ContiComp == 0L && length(branches) > 0L) {
    best        <- branches[which.max(CompSize[branches])]
    Class[best] <- 0L
  }
  
  Class
}