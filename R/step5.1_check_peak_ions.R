# check which MS1 ions are MS1 peaks from the peak list
#' Title
#'
#' @param spectrum.df spectrum
#' @param peaks.mat peaks matrix
#' @param rt_prec precursor rt
#'
#' @return spectrum with more information of the ions which are peaks also
#' @export
#'
check_ms1_ions <- function(spectrum.df, peaks.mat, rt_prec, rt_tol = 0.1){

  ions_peaks <- c()

  for (i in 1:nrow(spectrum.df)) {
    # i <- 60
    idx <- which( spectrum.df[i, 'mz_value'] >= peaks.mat[, 'mzmin'] &
                    spectrum.df[i, 'mz_value'] <= peaks.mat[, 'mzmax'] &
                    abs(rt_prec - peaks.mat[, 'rt']) < rt_tol )

    # print(paste(i, idx))
    if( length(unname(idx)) > 0 ){
      ions_peaks <- c(ions_peaks,  idx)
    }
  }

  return(ions_peaks)
}
