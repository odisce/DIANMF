# this script is to apply nGMCAs algorithm

#' update the matrix A (it contains elution profiles)
#'
#' @param Y matrix
#' @param A_init matrix
#' @param S matrix
#' @param lambda numeric sparsity constraint
#' @param options list of some parameters
#'
#' @return new_A matrix updated elution profiles
#'
updateA.f <- function(Y, A_init, S, lambda, options) {

  L <- svd( t(S) %*% S)$d[1]
  t <- 1
  R <- A_init
  A <- A_init

  for (i in 1:options$maxFBIteration) {
    # print(i)

    Anext <- pmax( R - (1/L) * (R %*% S - Y) %*% t(S), 0 )

    tnext <- ( 1+sqrt(1 + 4*(t^2)) ) /2
    R <- Anext + ((t-1)/tnext) * (Anext - A)

    if ( norm( Anext - A, 'F') / norm(A, 'F') < options$toleranceFB) {
      A <- Anext
      break }

    A <- Anext
    # print( error_function(Y, Anext, S) )
  }
  return(A)
}

#-------------------------------------------------------------------------------

#' apply soft-thresholding on a matrix
#'
#' @param X matrix
#' @param threshold numeric thresholding parameter
#'
#' @return matrix
#'
soft_threshold.f <- function(X, threshold){
  mat <- X
  mat <- abs(mat) - threshold
  mat <- pmax(mat, 0)
  mat <- sign(X) * mat
  return(mat)
}

#-------------------------------------------------------------------------------

#' update the matrix S (it contains pure spectra)
#'
#' @param Y matrix
#' @param A matrix
#' @param S_init matrix
#' @param lambda numeric sparsity constraint
#' @param options list of some parameters
#'
#' @return new_S matrix updated pure spectra
#'
updateS.f <- function(Y, A, S_init, lambda, options){

  L <- svd(A %*% t(A))$d[1]
  t <- 1
  R <- S_init
  S <- S_init

  for (i in 1:options$maxFBIteration) {
    # print(i)
    # lambda <- 0.9*lambda

    Snext <- pmax( soft_threshold.f( R - (1/L) * t(A) %*% (A %*% R - Y) , lambda/L) , 0)

    tnext <- ( 1+sqrt(1 + 4*(t^2)) ) /2

    R <- Snext + ((t-1)/tnext) * (Snext - S)

    if ( norm( Snext - S, 'F') / norm(S, 'F') < options$toleranceFB) {
      S <- Snext
      break }

    S <- Snext
    # print( error_function(Y, A, Snext) )
  }
  return(S)
}

#-------------------------------------------------------------------------------

#' calculate the difference between two matrices
#'
#' @param Y matrix
#' @param A matrix elution profiles
#' @param S matrix pure spectra
#'
#' @return numeric the error between Y and its approximation AS
#'
error_function <- function(Y, A, S){
  0.5*norm(Y - A %*% S, '2')^2 + norm(S, '1')
}

#-------------------------------------------------------------------------------

#'
#' apply nGMCAs on matrix Y to find A and S
#'
#' @param originalMatrix description
#' @param rank numeric rank of factorization = number of pure compounds in the mixed data
#' @param options list of parameters
#'
#' @return data: list of 2 matrices A and S
#'
nGMCAs <- function(originalMatrix = X_final, rank, options = options.l){

  if(options$useTranspose) { originalMatrix <- t(originalMatrix) }

  data <- list()
  res_nndsvd <- nndsvd(originalMatrix, rank)
  data$A <- res_nndsvd$W
  data$S <- res_nndsvd$H
  lambda <- 0.8 * max(originalMatrix)

  for (i in 1:options$maximumIteration) {
    # i <- 1
    print(i)
    data$S <- updateS.f(originalMatrix, data$A, data$S, lambda = lambda, options)

    data$A <- updateA.f(originalMatrix, data$A, data$S, lambda = lambda, options)

    print( error_function(originalMatrix, data$A, data$S ) )
  }

  if(options$useTranspose) {
    data$A <- t(data$A)
    data$S <- t(data$S)
  }
  return(data)
}

#-------------------------------------------------------------------------------
