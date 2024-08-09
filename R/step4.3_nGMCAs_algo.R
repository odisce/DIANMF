#' calculate the difference between Y an its approximation A %*% S
#'
#' @param Y matrix
#' @param A matrix elution profiles
#' @param S matrix pure spectra
#'
#' @return relative error
#' @export
error_function <- function(Y, A, S){
  norm(Y - A %*% S, '2') / norm(Y, '2')
}


#' update the matrix A (it contains elution profiles)
#'
#' @param Y matrix
#' @param A_init matrix
#' @param S matrix
#' @param maxFBIteration numeric max number of iteration
#' @param toleranceFB numeric stopping value
#'
#' @return new_A matrix updated elution profiles
#' @export
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


#' apply soft-thresholding on a matrix
#'
#' @param X matrix
#' @param threshold numeric thresholding parameter
#'
#' @return matrix
soft_threshold.f <- function(X, threshold){
  
  mat <- X
  mat <- abs(mat) - threshold
  mat <- pmax(mat, 0)
  mat <- sign(X) * mat
  
  return(mat)
}


#'  update the matrix S (it contains pure spectra)
#'
#' @param Y matrix
#' @param A matrix elution profiles
#' @param S_init matrix pure spectra
#' @param lambda numeric sparsity parameter of the pure spectra
#' @param maxFBIteration numeric FBS max number of iteration
#' @param toleranceFB numeric stopping value
#'
#' @return S
#' @export
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
      stop("NA values detected in Snext or S")
    }
    # Check for zero norm to avoid division by zero
    norm_S <- norm(S, 'F')
    if (norm_S == 0) {
      stop("Norm of S is zero, leading to division by zero")
    }
    if (norm(Snext - S, 'F') / norm_S < toleranceFB) {
      S <- Snext
      break
    }
    S <- Snext
  }
  
  return(S)
}


#' apply nGMCAs on matrix Y to find A and S
#'
#' @param originalMatrix description
#' @param rank numeric rank of factorization = number of pure compounds in the mixed data
#' @param errors_print Boolean 
#' @param initialization_method random, nndsvd or subSample 
#' @param maximumIteration numeric max number of iteration
#' @param maxFBIteration numeric FBS max number of iterations
#' @param toleranceFB numeric stopping value
#' @param H_sub matrix elution profiles
#'
#' @return data: list of 2 matrices A and S
#' @export
nGMCAs <- function(originalMatrix, rank,
                   maximumIteration = 10, maxFBIteration = 10, toleranceFB = 1e-5,
                  initialization_method = c('nndsvd', 'random', 'subSample'), H_sub = NULL,
                  errors_print = FALSE){

  originalMatrix <- base::t(originalMatrix)
  data <- list()
  
  if( initialization_method == 'nndsvd' ){
  res_nndsvd <- nndsvd_init( X = originalMatrix, rank = rank)
  data$A <- res_nndsvd$W
  data$S <- res_nndsvd$H
  } else if( initialization_method == 'random' ) {  
    res_random <- random_init(X = originalMatrix, rank = rank)
    data$A <- res_random$W
    data$S <- res_random$H
  } else if( initialization_method == 'subSample' & !is.null(H_sub) ){
    res_subSample <- subsample_init(Y = originalMatrix, rank = rank, H_sub)
    data$A <- res_subSample$A
    data$S <- res_subSample$S
  } else {
    print('error! H_sub must be !NULL when using the subSample initialization')
  }
  
  lambda <- 0.8 * max(originalMatrix)
  for (i in 1:maximumIteration) {
    
    data$S <- t(apply( data$S, 1, function(x) x / max(x)))
    data$S <- updateS.f( Y = originalMatrix, A = data$A, S_init = data$S, lambda = lambda,
                         maxFBIteration = maxFBIteration, toleranceFB = toleranceFB ) 

    data$A <- updateA.f(Y = originalMatrix, A_init = data$A, S = data$S,
                        maxFBIteration = maxFBIteration, toleranceFB = toleranceFB)
    
    if(errors_print) {print( error_function(originalMatrix, data$A, data$S )) }
  }
  
  data$A <- t(data$A)
  data$S <- t(data$S)

  return(data)
}
