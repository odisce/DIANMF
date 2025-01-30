# The MS2 pure spectrum that I got after NMF of the mixed MS2 data contains fragments that are related to:
# the precursor (peak) itself and its isotopes
# its adducts, losses, and their isotopes
# to other precursors of the same mz :):(

# Because of this, the pure MS2 spectrum should go to the second stage of filtering 
# i.e., to keep just the ions that exist in the isolation window where the precursor was fragmented


#' Get the total number of rows of previous matrices of a specific level in a list.
#'
#' @param my_list `list` of matrices.
#' @param level_index `numeric` index of the list level.
#'
#' @return `integer`
total_rows_before_level <- function(my_list, level_index) {
  if(level_index == 1){
    return(0)
  }
  
  if (level_index < 1 || level_index > length(my_list)) {
    stop("level_index is out of bounds")
  }
  
  total_rows <- 0
  for (i in 1:(level_index - 1)) {
    current_matrix <- my_list[[i]]
    if (is.matrix(current_matrix)) {
      total_rows <- total_rows + nrow(current_matrix)
    }
  }
  
  return(total_rows)
}


#' Extract the precursor fragments without their adducts and losses fragments.
#'
#' @param ms2_pure_spectrum `data.frame` MS2 pure spectrum related to the precursor.
#' @param ms2_matrices `list` of MS2 matrices.
#' @param mz_prec `numeric(1)` precursor mz.
#' @inheritParams extract_ms2_matrices
#'
#' @return `data.frame` MS2 spectrum of the precursor.
#' 
#' @export
filter_ms2_spectrum <- function(ms2_pure_spectrum, ms2_matrices, mz_prec, info.swath){

  idx.swath <- which(info.swath$lowerMz <= mz_prec & info.swath$upperMz >= mz_prec)
  if(length(idx.swath) > 1){
    print("warning, this peak belong to more than isolation window!")
    idx.swath <- idx.swath[1]   # it should be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!
  }
  
  if( length(idx.swath) > 0 ){
    start <- total_rows_before_level(my_list = ms2_matrices, level_index = idx.swath)
    end <- start + nrow(ms2_matrices[[idx.swath]])
    ms2_pure_spectrum_new <- ms2_pure_spectrum[ (start+1):end, ]
    # ms2_pure_spectrum_new$intensity <-  ms2_pure_spectrum_new$intensity / max( ms2_pure_spectrum_new$intensity)
    ms2_pure_spectrum_new <- ms2_pure_spectrum_new[ms2_pure_spectrum_new['intensity'] != 0, ]
    
    ms2_pure_spectrum_new <- as.data.frame(ms2_pure_spectrum_new)
    
  } else{  # if the mz_prec not included in any SWATH window, don't filter (this peaks wasn't fragmented in the MS2 level)
      print('No specific MS2 fragments.')
    return(NULL)
  }
 
  return(ms2_pure_spectrum_new)
}


# filter_ms2_spectrum <- function(ms2_pure_spectrum, ms2_matrices, mz_prec, info.swath, peak.idx = NULL){
#   
#   idx.swath <- which(info.swath$lowerMz <= mz_prec & info.swath$upperMz >= mz_prec)
#   if(length(idx.swath) > 1){
#     print("warning, this peak belong to more than isolation window!")
#     idx.swath <- idx.swath[1]
#   }
#   
#   if( length(idx.swath) > 0 ){
#     start <- total_rows_before_level(my_list = ms2_matrices, level_index = idx.swath)
#     end <- start + nrow(ms2_matrices[[idx.swath]])
#     ms2_pure_spectrum_new <- ms2_pure_spectrum
#     ms2_pure_spectrum_new[, "frag"] <- "na"
#     ms2_pure_spectrum_new[ (start+1):end, "frag"] <- "prec"
#     ms2_pure_spectrum_new$intensity <-  ms2_pure_spectrum_new$intensity / max( ms2_pure_spectrum_new$intensity)
#     ms2_pure_spectrum_new <- ms2_pure_spectrum_new[ms2_pure_spectrum_new['intensity'] != 0, ]
#     
#     ms2_pure_spectrum_new <- as.data.frame(ms2_pure_spectrum_new)
#     
#   } else{  # if the mz_prec not included in any SWATH window, don't filter (this peaks wasn't fragmented in the MS2 level)
#     if( !is.null(peak.idx) ){
#       print(paste(peak.idx, 'No specific MS2 fragments.'))
#     } else {
#       print('No specific MS2 fragments.')
#     }
#     return(NULL)
#   }
#   
#   return(ms2_pure_spectrum_new)
# }
