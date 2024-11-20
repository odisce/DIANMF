#' Calculate the difference between Y an its approximation A*S.
#'
#' @param Y `matrix`
#' @param A `matrix` every column is an elution profile.
#' @param S `matrix` every row is a pure spectrum.
#'
#' @return` numeric(1)` relative error.
error_function <- function(Y, A, S){
  norm(Y - A %*% S, '2') / norm(Y, '2')
}


#' Update A.
#'
#' @inheritParams error_function
#' @param A_init `matrix` every column is an elution profile.
#' @param maxFBIteration `numeric(1)` maximum number of iterations in the FB algorithm.
#' @param toleranceFB `numeric(1)` tolerance or stopping value of the FB algorithm.
#'
#' @return `matrix` A of updated elution profiles.
updateA.f <- function(Y, A_init, S, maxFBIteration, toleranceFB) {

  L <- svd( t(S) %*% S )$d[1]
  t <- 1
  R <- A_init
  A <- A_init

  for (i in 1:maxFBIteration) {

    Anext <- pmax( R - (1/L) * (R %*% S - Y) %*% t(S), 0 )
    tnext <- ( 1+sqrt(1 + 4*(t^2)) ) /2
    R <- Anext + ((t-1)/tnext) * (Anext - A)
    if ( norm( Anext - A, 'F') / norm(A, 'F') < toleranceFB) {
      A <- Anext
      break
    }
    A <- Anext
    # print( error_function(Y, Anext, S) )
  }
  
  return(A)
}


#' Apply soft-thresholding on a matrix.
#'
#' @inheritParams error_function
#' @param threshold `numeric(1)` thresholding parameter.
#'
#' @return `matrix`
soft_threshold.f <- function(Y, threshold){
  
  mat <- Y
  mat <- abs(mat) - threshold
  mat <- pmax(mat, 0)
  mat <- sign(Y) * mat
  
  return(mat)
}


#' Update S.
#'
#' @inheritParams error_function
#' @param S_init `matrix` every row is a pure spectrum.
#' @param lambda `numeric(1)` sparsity parameter of the pure spectra.
#' @inheritParams updateA.f
#'
#' @return `matrix` S of updated spectra.
updateS.f <- function(Y, A, S_init, lambda, maxFBIteration, toleranceFB){
 
  L <- svd(A %*% t(A))$d[1]
  t <- 1
  R <- S_init
  S <- S_init
  
  for (i in 1:maxFBIteration) {
    
    Snext <- pmax(soft_threshold.f(R - (1/L) * t(A) %*% (A %*% R - Y), lambda/L), 0)
    tnext <- (1 + sqrt(1 + 4 * (t^2))) / 2
    R <- Snext + ((t - 1) / tnext) * (Snext - S)
    # Check for NA values
    if (any(is.na(Snext)) || any(is.na(S))) {
      print("NA values detected in Snext or S")
      S <- Snext <- NULL
      break
    }
    # Check for zero norm to avoid division by zero
    norm_S <- norm(S, 'F')
    if (norm_S == 0) {
      print("Norm of S is zero, leading to division by zero")
      S <- Snext <- NULL
      break
    }
    if (norm(Snext - S, 'F') / norm_S < toleranceFB) {
      S <- Snext
      break
    }
    S <- Snext
  }
  
  return(S)
}


#' Approximate Y by A and S.
#'
#' @param X.m `matrix` contains components' sources mixed up in an unknown but linear way. Every row is an elution profile from the candidate precursor and its fragment ions through the rt scans in the columns.
#' @inheritParams random_init
#' @param errors_print `Logical` `TRUE` to print the error difference, `FALSE` otherwise. 
#' @param initialization_method `character` to specify the initialization method: random, nndsvd or subSample. 
#' @param maximumIteration `numeric(1)` maximum number of iterations.
#' @inheritParams updateA.f
#' @param H_sub `matrix` of elution profiles used with subSample initialization method.
#'
#' @return `list` of 2 matrices A and S.
#' 
#' @export
#' 
#' @examples
#'  m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3)
#'  
#' nGMCAs(X.m = m, rank = 2, maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
#'  initialization_method = 'nndsvd')
nGMCAs <- function(X.m, rank,
                   maximumIteration = 50, maxFBIteration = 100, toleranceFB = 1e-5,
                  initialization_method = c('nndsvd', 'random', 'subSample'), H_sub = NULL,
                  errors_print = FALSE){

  X.m <- base::t(X.m)
  data <- list()
  
  if( initialization_method == 'nndsvd' ){
  res_nndsvd <- nndsvd_init( X = X.m, rank = rank)
  data$A <- res_nndsvd$W
  data$S <- res_nndsvd$H
  } else if( initialization_method == 'random' ) {  
    res_random <- random_init(X = X.m, rank = rank)
    data$A <- res_random$W
    data$S <- res_random$H
  } else if( initialization_method == 'subSample' & !is.null(H_sub) ){
    res_subSample <- subsample_init(Y = X.m, rank = rank, H_sub)
    if( is.null(res_subSample) ){ return(NULL) }
    data$A <- res_subSample$A
    data$S <- res_subSample$S
  } else {
    print('error! H_sub must be !NULL when using the subSample initialization')
  }
  
  if( is.null(data$S) ){ return(NULL) }
  
  lambda <- 0.8 * max(X.m)
  for (i in 1:maximumIteration) {
    
    data$S <- t(apply( data$S, 1, function(x) x / max(x)))
    data$S <- updateS.f( Y = X.m, A = data$A, S_init = data$S, lambda = lambda,
                         maxFBIteration = maxFBIteration, toleranceFB = toleranceFB ) 
    
    data$A <- updateA.f(Y = X.m, A_init = data$A, S = data$S,
                        maxFBIteration = maxFBIteration, toleranceFB = toleranceFB)
    
    if(errors_print) {print( error_function(X.m, data$A, data$S )) }
  }
  
  data$A <- t(data$A)
  data$S <- t(data$S)

  return(data)
}
