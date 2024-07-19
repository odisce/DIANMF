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
random_init <- function(A, rank){
  
  max_val <- max(A)
  n_rows <- dim(A)[1]
  n_cols <- dim(A)[2]
  
  W <- matrix( runif(n_rows * rank, min = 0, max = max_val), nrow = n_rows, ncol = rank) 
  H <- matrix( runif(rank * n_cols, min = 0, max = max_val), nrow = rank, ncol = n_cols) 
    
  return(list(
    'W' = W,
    'H' = H
  ))
}
#-------------------------------------------------------------------------------

## sub-sampling
# to come soon
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
nndsvd <- function(A, k){
  
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
