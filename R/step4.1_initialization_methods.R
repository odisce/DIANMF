## random initialization

#' Initialize matrices W and H randomly.
#
#' @param X matrix
#' @param rank of factorization
#'
#' @return W and H random matrices
#' @export
#' @importFrom stats runif
random_init <- function(X, rank){
  
  max_val <- max(X)
  n_rows <- dim(X)[1]
  n_cols <- dim(X)[2]
  W <- matrix( runif(n_rows * rank, min = 0, max = max_val), nrow = n_rows, ncol = rank) 
  H <- matrix( runif(rank * n_cols, min = 0, max = max_val), nrow = rank, ncol = n_cols) 
    
  return(list(
    'W' = W,
    'H' = H ))
}
#-------------------------------------------------------------------------------

##  subSample initialization: given the elution profiles, initialize S

#' subSample initialization
#'
#' @param Y matrix
#' @param rank of factorization
#' @param H_sub pure MS1 elution profiles matrix
#'
#' @return A and S initializing matrices based on the data
#' @export
subsample_init <- function(Y, rank, H_sub){  # This method is just for NMF of MS2 data. (and can't be applied for MS1 data)
  
  # for S, solve the inverse problem: argmin_{S >=0} 0.5||Y - AS||_F^2 + lambda*||S||_1 (by FBS).
  # max_val <- max(Y) 
  S <- matrix( runif( rank * ncol(Y), min = 0, max = 1), nrow = rank, ncol = ncol(Y)) # initialize S
  
  lambda <- 0.8 * max(Y)
  S <- updateS.f( Y, A = t(H_sub), S_init = S, lambda = lambda,
                       maxFBIteration = 10, toleranceFB = 1e-5 ) 
  
  return(list(
    'A' = t(H_sub),
    'S' = S ))
}
#-------------------------------------------------------------------------------

## nndsvd method non-negative double singular value decomposition

#' positive projection
#'
#' @param x vector of numerics
#'
#' @return vector the projection of x in the positive orthant
pos <- function(x){ as.numeric(x>=0) * x }


#' negative projection
#'
#' @param x vector of numerics
#'
#' @return vector the projection of x in the negative orthant
neg <- function(x){ as.numeric(x<=0) * abs(x) }


#' calculates the square root of the sum of the squared elements in the input matrix "x"
#'
#' @param x matrix
#'
#' @return numeric
norm_euc <- function(x){ sqrt(drop(crossprod(x))) }


#' apply nndsvd (non-negative double singular value decomposition) on a matrix
#'
#' @param X matrix
#' @param rank of factorization
#'
#' @return list of two matrices H and W
#' @export
nndsvd_init <- function(X, rank){
  
  if( any(X<0) ){
    print('The input matrix contains negative elements !')
    return(NULL)
  }
  
  #the matrices of the factorization
  W = matrix(0, dim(X)[1], rank);
  H = matrix(0, rank, dim(X)[2]);
  
  svd_res = svd(X, rank, rank);
  U <- svd_res$u; S <- svd_res$d; V <- svd_res$v
  
  nb <- min(rank, dim(V)[2])
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
    "H" = H ))
}
#-------------------------------------------------------------------------------
