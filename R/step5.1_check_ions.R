#' Check which ions are also detected MS1 peak.
#'
#' @param ions_maybe_peaks.v `numeric` vector, ions mz values.
#' @inheritParams extract_ms_matrix.f
#' @param rt_prec `numeric(1)` apex retention time.
#' @param rt_tol `numeric(1)` retention time tolerance.
#'
#' @return `numeric` MS1 peaks indexes.
check_ms1_ions <- function(ions_maybe_peaks.v, ms1_peaks.df, rt_prec, rt_tol = 2){
  ions_peaks <- c()
  for (i in 1:length(ions_maybe_peaks.v)) {
    idx <- which( ions_maybe_peaks.v[i] >= ms1_peaks.df[, 'mzmin'] & 
                    ions_maybe_peaks.v[i] <= ms1_peaks.df[, 'mzmax'] &
                    abs(rt_prec - ms1_peaks.df[, 'rt']) < rt_tol )
    if( length(unname(idx)) > 0 ){
      ions_peaks <- c(ions_peaks,  idx)
    }
  }
  
  return(ions_peaks)
}
