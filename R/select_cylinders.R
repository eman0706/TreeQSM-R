#' Select cylinders by index
#'
#' @param cylinder Cylinder structure
#' @param Ind Index of cylinders to keep
#' @return Updated cylinder structure
select_cylinders <- function(cylinder, Ind) {
  Names <- names(cylinder)
  n <- length(Names)
  for (i in 1:n) {
    if (is.matrix(cylinder[[Names[i]]])) {
      cylinder[[Names[i]]] <- cylinder[[Names[i]]][Ind, , drop = FALSE]
    } else if (is.vector(cylinder[[Names[i]]]) || is.factor(cylinder[[Names[i]]])) {
      cylinder[[Names[i]]] <- cylinder[[Names[i]]][Ind]
    }
  }
  return(cylinder)
}
