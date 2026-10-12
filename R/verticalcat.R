#' Verticalcat
#'
#' Equivalente di verticalcat.m: concatena i cover set di tutti i layer
#' di un segmento in un unico vettore "Sets", registrando in "IndSets"
#' (nl x 2) l'indice di inizio/fine di ciascun layer dentro Sets.
#' @param Seg Lista di layer (ciascuno un vettore di indici di cover set)
#' @return Lista con Sets (vettore concatenato) e IndSets (nl x 2)
#' @export
verticalcat <- function(Seg) {
  nl <- length(Seg)
  lens <- sapply(Seg, length)
  ends <- cumsum(lens)
  starts <- c(1L, if (nl > 1) head(ends, -1) + 1L else integer(0))
  Sets <- unlist(Seg)
  IndSets <- cbind(starts, ends)
  list(Sets = Sets, IndSets = IndSets)
}
