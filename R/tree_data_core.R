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

#' Tree Data Core
#'
#' Porting fedele PARZIALE di tree_data.m (Version 3.0.1, Copyright (C)
#' 2013-2022 Pasi Raumonen). Include solo le parti direttamente rilevanti
#' al problema affrontato in sessione (identificazione del tronco,
#' volume, altezza, DBH) - NON include crown_measures, triangulate_stem,
#' ne' le numerose distribuzioni (per diametro/altezza/angolo/ordine di
#' ramo), che sono fuori dallo scopo attuale.
#'
#' NOTA IMPORTANTE (scoperta durante il porting): il vero algoritmo
#' identifica il tronco con "Trunk = cylinder.branch == 1" preso in
#' ORDINE DI INDICE (non con una ricostruzione di catena come la nostra
#' identify_stem_chain() - quella era una nostra invenzione per
#' compensare la mancanza di questo file). Questo significa che il vero
#' algoritmo si aspetta che TUTTI i cilindri del tronco fisico, in
#' ordine di creazione, abbiano branch==1 - se questo non vale sui
#' nostri dati (tronco frammentato su piu' numeri di branch), il
#' problema e' a monte (probabile residuo in segments.m/cylinders.m),
#' non qualcosa da aggirare qui.
#'
#' Il DBH vero ("DBHqsm") usa la LUNGHEZZA CUMULATIVA lungo l'asse dei
#' cilindri del tronco (non l'altezza verticale, come avevamo assunto in
#' un fix precedente della sessione) - fedele all'originale. Il vero
#' algoritmo poi rifitta un cilindro DEDICATO sui punti reali della
#' nuvola in una fascia stretta 1.1-1.5m, con verifica di plausibilita',
#' invece di fidarsi ciecamente del cilindro grezzo trovato a 1.3m -
#' questo e' il meccanismo reale per gestire basi rumorose/contaminate,
#' molto piu' robusto delle euristiche di correzione a valle inventate
#' in sessione.
#'
#' @param cylinder Output di cylinders() (radius, length, start, axis, branch)
#' @param trunk_points Nuvola di punti REALI del tronco (matrice Nx3) -
#'   serve per il rifit dedicato del DBH
#' @return Lista con i campi principali di treedata (TotalVolume,
#'   TrunkVolume, BranchVolume, TreeHeight, TrunkLength, BranchLength,
#'   TotalLength, TrunkArea, BranchArea, TotalArea, DBHqsm, DBHcyl,
#'   Trunk = indici dei cilindri di tronco in ordine di indice)
#' @export
tree_data_core <- function(cylinder, trunk_points = NULL) {
  Rad <- cylinder$radius
  Len <- cylinder$length
  nc  <- length(Rad)
  ind <- seq_len(nc)
  Trunk <- cylinder$branch == 1L
  
  treedata <- list()
  treedata$TotalVolume  <- 1000 * pi * sum(Rad^2 * Len)
  treedata$TrunkVolume  <- 1000 * pi * sum(Rad[Trunk]^2 * Len[Trunk])
  treedata$BranchVolume <- 1000 * pi * sum(Rad[!Trunk]^2 * Len[!Trunk])
  
  bottom <- min(cylinder$start[, 3])
  i_top  <- which.max(cylinder$start[, 3])
  top    <- cylinder$start[i_top, 3]
  if (cylinder$axis[i_top, 3] > 0) {
    top <- top + Len[i_top] * cylinder$axis[i_top, 3]
  }
  treedata$TreeHeight <- top - bottom
  
  treedata$TrunkLength  <- sum(Len[Trunk])
  treedata$BranchLength <- sum(Len[!Trunk])
  treedata$TotalLength  <- treedata$TrunkLength + treedata$BranchLength
  
  treedata$TrunkArea  <- 2 * pi * sum(Rad[Trunk] * Len[Trunk])
  treedata$BranchArea <- 2 * pi * sum(Rad[!Trunk] * Len[!Trunk])
  treedata$TotalArea  <- 2 * pi * sum(Rad * Len)
  
  dbh <- dbh_cylinder(Trunk, cylinder, ind, trunk_points)
  treedata$DBHqsm <- dbh$DBHqsm
  treedata$DBHcyl <- dbh$DBHcyl
  treedata$Trunk  <- ind[Trunk]
  
  treedata
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

#' DBH Cylinder
#'
#' Porting fedele di dbh_cylinder (sotto-funzione di tree_data.m).
#' Calcola il DBH grezzo dal QSM (DBHqsm, lunghezza cumulativa lungo
#' l'asse dei cilindri di tronco) e, se disponibili abbastanza punti
#' reali nella fascia 1.1-1.5m, rifitta un cilindro dedicato (DBHcyl),
#' accettandolo solo se supera i controlli di plausibilita'
#' dell'originale.
#'
#' @param Trunk Vettore logico (branch==1), allineato a cylinder
#' @param cylinder Output di cylinders()
#' @param ind Indici 1:nc
#' @param trunk_points Nuvola di punti reali del tronco (matrice Nx3),
#'   o NULL se non disponibile (in tal caso DBHcyl = DBHqsm, fedele al
#'   ramo "else" dell'originale per pochi punti)
#' @export
dbh_cylinder <- function(Trunk, cylinder, ind, trunk_points = NULL) {
  n <- sum(Trunk)
  T <- ind[Trunk]
  if (n == 0L || length(T) == 0L) {
    return(list(DBHqsm = NA_real_, DBHcyl = NA_real_))
  }
  
  i <- 1L
  while (i < n && sum(cylinder$length[T[seq_len(i)]]) < 1.3) {
    i <- i + 1L
  }
  DBHqsm <- 2 * cylinder$radius[T[i]]
  
  DBHcyl <- DBHqsm
  if (!is.null(trunk_points) && nrow(trunk_points) > 0) {
    V <- sweep(trunk_points, 2, cylinder$start[1, ])
    h <- as.vector(V %*% cylinder$axis[1, ])
    I <- h < 1.5 & h > 1.1
    if (sum(I) > 100L) {
      Tpts <- trunk_points[I, , drop = FALSE]
      # NOTA FEDELTA': l'originale usa "i" (la posizione all'interno del
      # sottoinsieme di tronco), non "T[i]" (l'indice reale nell'array
      # completo dei cilindri) - sia qui che nel controllo di
      # allineamento sotto. Ha senso solo se i cilindri di tronco
      # occupano esattamente le prime n posizioni dell'array (il
      # segmento 1 cattura da solo l'intero tronco fisico, scritto per
      # primo e in blocco) - ulteriore conferma, insieme a "Trunk =
      # cylinder.branch==1" preso in ordine di indice, che il vero
      # algoritmo si aspetta il tronco intero dentro il segmento 1.
      # Preservato cosi' com'e', non "corretto" a T[i].
      cyl0 <- select_cylinders(cylinder, i)
      cyl  <- tryCatch(least_squares_cylinder(Tpts, cyl0), error = function(e) NULL)
      if (!is.null(cyl)) {
        RadiusOK <- (2 * cyl$radius > 0.8 * DBHqsm) && (2 * cyl$radius < 1.2 * DBHqsm)
        axis_align <- abs(sum(cylinder$axis[i, ] * cyl$axis))
        if (isTRUE(RadiusOK) && is.finite(axis_align) && axis_align > 0.9 &&
            isTRUE(cyl$conv) && isTRUE(cyl$rel)) {
          DBHcyl <- 2 * cyl$radius
        }
      }
    }
  }
  
  list(DBHqsm = DBHqsm, DBHcyl = DBHcyl)
}
