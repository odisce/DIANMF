#' Given X.m and W, find H.
#'
#' @inheritParams nGMCAs
#' @param W pure spectra `matrix`.
#'
#' @return H: elution profiles (chromatograms) `matrix`.
#' @export
find_H <- function(X.m, W, maxFBIteration = 20, toleranceFB = 1e-5){
  
  A <- matrix( runif(ncol(X.m) * ncol(W), min = 0, max = 1), nrow = ncol(W), ncol = ncol(X.m))
  A <- updateA.f(Y = t(X.m), A_init = t(A), S = t(W), maxFBIteration = maxFBIteration, toleranceFB = toleranceFB)
  
  return(t(A))
}


#' Given X.m and H, find W.
#'
#' @inheritParams nGMCAs
#' @param H elution profiles `matrix`.
#'
#' @return W: pure spectra `matrix`.
#' @export
find_W <- function(X.m, H, maxFBIteration = 20, toleranceFB = 1e-5){
  
  W <- matrix( runif( nrow(X.m) * nrow(H), min = 0, max = 1), nrow = nrow(X.m), ncol = nrow(H))
  lambda <- 0.8 * max(X.m)
  W <- updateS.f(Y = t(X.m), A = t(H), S_init = t(W), lambda = lambda, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB)
  
  # normalize every row in W by the the row's max value
  row_max <- apply(W, 1, max)
  W <- W / row_max
  
  return(t(W))
}
