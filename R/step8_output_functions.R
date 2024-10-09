#' Extract pure spectra.
#'
#' @param features feature list. 
#' @param spec_level `character` determine the pure spectrum spec_level 'MS1', 'MS2' or 'MS2_specific'. 'MS2_specific' is the spectrum related to the precursor and its fragments without its adducts and loss fragments.
#'
#' @return `list` of pure spectra.
#' 
#' @export
extract_pureSpect <- function(features, spec_level){
  
  if( spec_level == "MS1" ){
    spect <- lapply(features, function(x) x$MS1_pure_spectrum)
  } else if( spec_level == "MS2" ){
    spect <- lapply(features, function(x) x$MS2_pure_spectrum)
    
  } else {
    spect <- lapply(features, function(x) x$MS2_pure_spectrum_specific)
  }
  
  if( length(spect) == 1){
    spect <- spect[[1]]
  }
  
  return(spect)
}


# extrcat_pureChrom <- function(features, ms_level, comp_idx){
#   
#   if( ms_level == "MS1" ){
#     chrom <- lapply(features, function(x) x$H_ms1)
#   } else {
#     chrom <- lapply(features, function(x) x$H_ms2)
#   }
#   
#   if( !is.null(comp_idx) ){
#     chrom <- lapply(chrom, function(x) x[comp_idx, ])
#   }
#   
#   return(chrom)
# }


#' Extract raw (mixed) matrices.
#'
#' @inheritParams extract_pureSpect 
#' @param ms_level `character` 'MS1' or 'MS2' level.
#'
#' @return `list` of raw matrices before factorization.
#' 
#' @export
extract_mixedMat <- function(features, ms_level) {
  
  if (ms_level == "MS1") {
    mat <- lapply(features, function(x) x$MS1_mixed_mat)
  } else {
    mat <- lapply(features, function(x) do.call(rbind,x$MS2_mixed_mat))
  }
  
  if( length(mat) == 1){
    mat <- mat[[1]]
  }
  
  return(mat)
}


#' Extract factorization matrices.
#'
#' @inheritParams extract_mixedMat
#' @param H `Logical` if H is TRUE extract H matrices else extract W.  
#'
#' @return `list` of H or W matrices.
#' 
#' @export
extract_pureMat <- function(features, ms_level, H){
  
  if (ms_level == "MS1") {
    if ( H == TRUE ){
      mat <- lapply(features, function(x) x$H_ms1)
    } else {
      mat <- lapply(features, function(x) x$W_ms1)
    }
  } else {
    if ( H == TRUE ){
      mat <- lapply(features, function(x) x$H_ms2)
    } else {
      mat <- lapply(features, function(x) x$W_ms2)
    }
  }
  
  if( length(mat) == 1){
    mat <- mat[[1]]
  }
  
  return(mat)
}
