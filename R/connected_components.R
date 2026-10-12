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

connected_components <- function(Nei, Sub, MinSize, Fal = NULL) {
  
  # ---------------------------------------------------------------------
  # CONNECTED_COMPONENTS.R      Determines the connected components of cover
  #                             sets using their neighbour-relation
  #
  # Version 1.0
  # Latest update     2 Mar 2026
  #
  # Copyright (C) 2013-2017 Pasi Raumonen
  # ---------------------------------------------------------------------
  
  # Determines connected components of the subset of cover sets defined
  # by "Sub" such that each component has at least "MinSize"
  # number of cover sets.
  #
  # Inputs:
  # Nei       Neighboring cover sets of each cover set, list of vectors
  # Sub       Subset whose components are determined,
  #               length(Sub) < 2 means no subset and thus the whole point cloud
  #               "Sub" may be also a vector of cover set indexes in the subset
  #               or a logical vector, where n_sets is the number of
  #               all cover sets
  # MinSize   Minimum number of cover sets in an acceptable component
  # Fal       Logical false vector for the cover sets
  #
  # Outputs:
  # Components    Connected components, list of integer vectors
  # CompSize      Number of sets in the components, integer vector
  
  # Handle small subsets (at most 3 cover sets)
  n <- length(Sub)
  
  if (n <= 3 && !is.null(Sub) && !any(is.na(Sub)) && length(Sub) > 0 && Sub[1] > 0) {
    if (n == 1) {
      Components <- list(as.integer(Sub))
      CompSize <- 1L
    } else if (n == 2) {
      I <- Nei[[Sub[1]]] == Sub[2]
      if (any(I)) {
        Components <- list(as.integer(Sub))
        CompSize <- 1L
      } else {
        Components <- list(as.integer(Sub[1]), as.integer(Sub[2]))
        CompSize <- c(1L, 1L)
      }
    } else if (n == 3) {
      I <- Nei[[Sub[1]]] == Sub[2]
      J <- Nei[[Sub[1]]] == Sub[3]
      K <- Nei[[Sub[2]]] == Sub[3]
      num_connections <- sum(c(any(I), any(J), any(K)))
      
      if (num_connections >= 2) {
        Components <- list(as.integer(Sub))
        CompSize <- 1L
      } else if (any(I)) {
        Components <- list(as.integer(Sub[1:2]), as.integer(Sub[3]))
        CompSize <- c(2L, 1L)
      } else if (any(J)) {
        Components <- list(as.integer(Sub[c(1, 3)]), as.integer(Sub[2]))
        CompSize <- c(2L, 1L)
      } else if (any(K)) {
        Components <- list(as.integer(Sub[2:3]), as.integer(Sub[1]))
        CompSize <- c(2L, 1L)
      } else {
        Components <- list(as.integer(Sub[1]), as.integer(Sub[2]), as.integer(Sub[3]))
        CompSize <- c(1L, 1L, 1L)
      }
    }
    
  } else if (any(Sub) || (length(Sub) == 1 && Sub[1] == 0)) {
    nb <- length(Nei)
    
    Fal_not_given <- is.null(Fal)
    if (Fal_not_given) {
      Fal <- rep(FALSE, nb)
    }
    
    if (length(Sub) == 1 && Sub[1] == 0) {
      # All the cover sets
      ns <- nb
      if (Fal_not_given) {
        Sub <- rep(TRUE, nb)
      } else {
        Sub <- !Fal
      }
    } else if (!is.logical(Sub)) {
      # Subset of cover sets
      ns <- length(Sub)
      if (Fal_not_given) {
        sub <- rep(FALSE, nb)
      } else {
        sub <- Fal
      }
      sub[Sub] <- TRUE
      Sub <- sub
    } else {
      # Subset of cover sets
      ns <- sum(Sub)
    }
    
    Components <- vector("list", ns)
    CompSize <- integer(ns)
    nc <- 0
    m <- which(Sub)[1]
    
    i <- 0
    Comp <- integer(ns)
    
    while (i < ns) {
      Add <- Nei[[m]]
      I <- Sub[Add]
      Add <- Add[I]
      a <- length(Add)
      Comp[1] <- m
      Sub[m] <- FALSE
      t <- 1
      
      while (a > 0) {
        Comp[(t + 1):(t + a)] <- Add
        Sub[Add] <- FALSE
        t <- t + a
        Add <- unlist(Nei[Add], use.names = FALSE)
        I <- Sub[Add]
        Add <- Add[I]
        
        # Select the unique elements of Add
        n_add <- length(Add)
        if (n_add > 2) {
          I <- rep(TRUE, n_add)
          for (j in 1:n_add) {
            if (Fal[Add[j]]) {
              I[j] <- FALSE
            } else {
              Fal[Add[j]] <- TRUE
            }
          }
          Fal[Add] <- FALSE
          Add <- Add[I]
        } else if (n_add == 2) {
          if (Add[1] == Add[2]) {
            Add <- Add[1]
          }
        }
        a <- length(Add)
      }
      
      i <- i + t
      if (t >= MinSize) {
        nc <- nc + 1
        Components[[nc]] <- as.integer(Comp[1:t])
        CompSize[nc] <- t
      }
      
      if (i < ns) {
        remaining <- which(Sub)
        if (length(remaining) > 0) {
          m <- remaining[1]
        } else {
          break
        }
      }
    }
    
    Components <- Components[1:nc]
    CompSize <- CompSize[1:nc]
    
  } else {
    Components <- list()
    CompSize <- integer(0)
  }
  
  return(list(Components = Components, CompSize = CompSize))
}