#' Average
#'
#' Equivalente di average.m: media per colonna di una matrice Nx3 (o di
#' un singolo punto, che viene restituito invariato).
#' @param P Matrice n x 3, o vettore lunghezza 3
#' @return Vettore lunghezza 3
#' @export
average <- function(P) {
  if (is.null(dim(P))) return(as.numeric(P))
  colMeans(P)
}
