#' check which MS1 ions are MS1 peaks from the peak list
#'
#' @param ions_maybe_peaks.v numeric_vector fragments mz values
#' @param peaks.df data_frame of all peaks
#' @param rt_prec numeric(1) peak/precursor retention time at the apex
#' @param rt_tol numeric(1) retention time tolerance
#'
#' @return vector of peaks which are also ions
check_ms1_ions <- function(ions_maybe_peaks.v, peaks.df, rt_prec, rt_tol = 2){
  ions_peaks <- c()
  for (i in 1:length(ions_maybe_peaks.v)) {
    idx <- which( ions_maybe_peaks.v[i] >= peaks.df[, 'mzmin'] & 
                    ions_maybe_peaks.v[i] <= peaks.df[, 'mzmax'] &
                    abs(rt_prec - peaks.df[, 'rt']) < rt_tol )
    if( length(unname(idx)) > 0 ){
      ions_peaks <- c(ions_peaks,  idx)
    }
  }
  
  return(ions_peaks)
}
