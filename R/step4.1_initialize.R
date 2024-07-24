# this script contains initialization methods

## random initialization


#' Initialize matrices W and H randomly.
#'
#' @param A 
#' @param rank 
#'
#' @return W and H randon matrices
#' @export
#' @importFrom stats runif
#'
random_init <- function(A, k){
  
  max_val <- max(A)
  n_rows <- dim(A)[1]
  n_cols <- dim(A)[2]
  
  W <- matrix( runif(n_rows * k, min = 0, max = max_val), nrow = n_rows, ncol = k) 
  H <- matrix( runif(k * n_cols, min = 0, max = max_val), nrow = k, ncol = n_cols) 
    
  return(list(
    'W' = W,
    'H' = H
  ))
}
#-------------------------------------------------------------------------------

## initialize the elution profiles matrix by the pure ones I got from NMF MS1 data

#' subSample initialization
#'
#' @param A 
#' @param k 
#' @param H_sub pure MS1 elution profiles
#'
#' @return W and H initializing matrices based on the data
#' @export
subsample_init <- function(mat, k, H_sub){  # This method is for NMF of MS2 data.
  
  # for W, solve the inverse problem: argmin_{W >=0} 0.5||A - WH||_F^2 (by FBS, specifically gradient projection algorithm).
  max_val <- max(mat)/20
  W <- matrix( runif( k * ncol(mat), min = 0, max = max_val), nrow = k, ncol = ncol(mat)) 
  param <- list(
    'maxFBIteration' = 10,
    'toleranceFB' = 1e-5
  )
  
  L <- svd(H_sub %*% t(H_sub))$d[1]
  t <- 1
  R <- W
  S <- W
  
  for (i in 1:20) {
    # print(i)
    Snext <- pmax( R - (1/L) * t(H_sub) %*% (H_sub %*% R - mat) , 0)
    
    tnext <- ( 1+sqrt(1 + 4*(t^2)) ) /2
    
    R <- Snext + ((t-1)/tnext) * (Snext - S)
    
    if ( norm( Snext - S, 'F') / norm(S, 'F') < 1e-5) {
      S <- Snext
      break }
    
    S <- Snext
    # print( error_function(mat, H_sub, Snext) )
  }
  
  return(list(
    'W' = S,
    'H' = H_sub
  ))
}
#-------------------------------------------------------------------------------

# nndsvd method non-negative double singular value decomposition

#' positive projection
#'
#' @param x vector of numeric
#'
#' @return vector the projection of x in the positive orthant
#'
pos <- function(x){ as.numeric(x>=0) * x }

#' negative projection
#'
#' @param x vector of numeric
#'
#' @return vector the projection of x in the negative orthant
#'
neg <- function(x){ as.numeric(x<=0) * abs(x) }

#' calculates the square root of the sum of the squared elements in the input matrix "x"
#'
#' @param x matrix
#'
#' @return numeric
norm_euc <- function(x){ sqrt(drop(crossprod(x))) }

# A <- mixed_eics
# k <- 3

#' apply nndsvd on a matrix
#'
#' @param A matrix
#' @param k integer: the rank
#'
#' @return list of two matrices H and W
#'
nndsvd_init <- function(A, k){
  
  if( any(A<0) ){
    print('The input matrix contains negative elements !')
    return(NULL)
  }
  
  #the matrices of the factorization
  W = matrix(0, dim(A)[1], k);
  H = matrix(0, k, dim(A)[2]);
  
  svd_res = svd(A, k, k);
  U <- svd_res$u; S <- svd_res$d; V <- svd_res$v
  
  nb <- min(k, dim(V)[2])
  for( i in 1:nb ){
    uu = U[,i]; vv = V[,i];
    
    uPositive <- pos(uu)
    vPositive <- pos(vv)
    
    uPositiveNorm <- norm_euc(uPositive)
    vPositiveNorm <- norm_euc(vPositive)
    mp <- uPositiveNorm * vPositiveNorm
    
    uNegative <- neg(uu)
    vNegative <- neg(vv)
    
    uNegativeNorm <- norm_euc(uNegative)
    vNegativeNorm <- norm_euc(vNegative)
    mn <- uNegativeNorm * vNegativeNorm
    
    if (mp >= mn){
      W[,i] = sqrt(S[i] * mp) * uPositive / uPositiveNorm
      H[i,] = sqrt(S[i] * mp) * vPositive / vPositiveNorm;
    }else{
      W[,i] = sqrt(S[i] * mn) * uNegative / uNegativeNorm;
      H[i,] = sqrt(S[i] * mn) * vNegative / vNegativeNorm;
    }
  }
  
  return(list(
    "W" = W,
    "H" = H))
}
#-------------------------------------------------------------------------------
