# PROVA: quanto conta la base ortonormale casuale dell'originale?
# L'originale (orthonormal_vectors.m) parte da un vettore casuale rand(3,1); il porting parte da un vettore di riferimento fisso. La base ortonormale fissa l'origine angolare
# della griglia di settori della copertura di superficie (surface_coverage). Questa prova sostituisce, nella sessione pulita, orthonormal_vectors con una versione che parte da un
# vettore pseudo-casuale diverso a ogni chiamata, SENZA toccare il generatore di numeri casuali di R (altrimenti si cambierebbero anche tutti gli altri numeri casuali), e
# riesegue gli alberi 4, 16 e 35 (seed 0) confrontando con le corse validate (variante di average: colMeans, come nelle corse validate).
# PREVISIONI SCRITTE PRIMA
#  P_O1: tutti i diametri (DBHcyl ufficiale e diametri alle tre quote) restano entro 0.05 cm in almeno il 95% dei 48 modelli: l'effetto della base fissa e' trascurabile sui diametri.
#  P_O2: i volumi fino a 4 m restano entro 0.1 litri in almeno il 90% dei modelli.
#  Se P_O1 non si conferma, i risultati dipendono dall'origine angolare dei settori (un elemento stocastico dell'originale che il porting non riproduce): e' un risultato da riportare.
# USO (sessione pulita, dopo carica_libreria_pulita()):  source("prova_orthonormal_casuale.R");  r <- prova_orthonormal_casuale()
prova_orthonormal_casuale <- function(base = "F:/EB/tartarino/sw python/TreeQSM_ProjectR/R", alberi = c(4, 16, 35), nuvole = c("alberi_corretti", "alberi_corretti_las11")) {
  stopifnot(exists("riproduci_pulita"), exists("orthonormal_vectors"), exists("cross_product"))
  orig <- get("orthonormal_vectors", envir = globalenv()); contatore <- 0L
  casuale <- function(U) {
    contatore <<- contatore + 1L; k <- contatore
    V0 <- c(sin(12.9898 * k), sin(78.233 * k + 1), sin(37.719 * k + 2))          # vettore pseudo-casuale, indipendente dal generatore di R
    V <- cross_product(V0, U); nV <- sqrt(sum(V^2))
    if (!is.finite(nV) || nV < 1e-10) { V <- cross_product(c(1, 0, 0), U); nV <- sqrt(sum(V^2)) }
    W <- cross_product(V, U); W <- W / sqrt(sum(W^2)); V <- V / nV
    list(V = V, W = W)
  }
  assign("orthonormal_vectors", casuale, envir = globalenv()); on.exit(assign("orthonormal_vectors", orig, envir = globalenv()), add = TRUE)
  cat("orthonormal_vectors sostituita da una versione con orientamento pseudo-casuale; riesecuzione di", paste(alberi, collapse = ", "), "\n")
  riproduci_pulita(base = base, alberi = alberi, nuvole = nuvole)
}
