#' Check which ions are also detected as MS1 peaks.
#'
#' @inheritParams extract_ms1_pure_spectrum
#' @param comp_ms1 index of the corresponding peak compound.
#' @inheritParams extract_ms_matrix.f
#' @param rt_prec `numeric(1)` apex retention time.
#' @param rt_tol `numeric(1)` retention time tolerance.
#'
#' @return `numeric` MS1 peaks indexes.
check_ms1_ions <- function(W_ms1, comp_ms1, ms1_peaks.df, rt_prec, rt_tol = 5){
  
  ions_maybe_peaks <- which( W_ms1[,comp_ms1] >= 0.7 * rowSums(W_ms1) )
  # ions_maybe_peaks <- which( W_ms1[,comp_ms1] >= rowSums(W_ms1[, -comp_ms1]) )
  ions_maybe_peaks <- as.numeric(names(ions_maybe_peaks))

  ions_peaks <- c()
  for (i in 1:length(ions_maybe_peaks)) {
    idx <- which( ions_maybe_peaks[i] >= ms1_peaks.df[, 'mzmin'] & 
                    ions_maybe_peaks[i] <= ms1_peaks.df[, 'mzmax'] &
                    abs(rt_prec - ms1_peaks.df[, 'rt']) < rt_tol )
    if( length(unname(idx)) > 0 ){
      ions_peaks <- c(ions_peaks,  idx)
    }
  }
  
  return(ions_peaks)
}
