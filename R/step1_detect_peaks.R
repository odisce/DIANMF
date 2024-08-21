#' Detect MS1 peaks using XCMS.
#'
#' @param rawData.onDiskMSnExp `OnDiskMSnExp` object for `onDisk` mode.
#' @param ... Additional parameters passed to \code{\link[xcms]{CentWaveParam}}.
#'
#' @return MS1 peaks `matrix`.
#' @export
#' @importFrom xcms findChromPeaks CentWaveParam chromPeaks
detect_peaks_by_xcms <- function(rawData.onDiskMSnExp, ...){
  
  # create centwave parameter object
  cwp <- xcms::CentWaveParam(...)
  # detect EIC peaks
  detected_eics <- xcms::findChromPeaks(rawData.onDiskMSnExp, param = cwp)
  
  ms1_peaks.mat <- xcms::chromPeaks(detected_eics)
  
  return(ms1_peaks.mat)
}


#' Prepare the MS1 peaks.
#' 
#' @param ms1_peaks `matrix` or `data.frame` contains MS1 peaks to be identified. It should contains these columns: mz, mzmin, mzmax, rt, rtmin, rtmax, into.
#'
#' @return MS1 peaks as `data.frame`, ordered in decreasing order of 'into' column, and column is_ion is added.
#' @export
#' @importFrom dplyr arrange desc
#' @import magrittr
prepare_ms1_peaks <- function(ms1_peaks){
  
  if( !("matrix" %in% class(ms1_peaks)) & !("data.frame" %in% class(ms1_peaks)) ){
    print("Error; the ms1_peaks should be a matrix or a data.frame R objects")
    return(NULL) }
  
  if( "matrix"  %in% class(ms1_peaks) ){
    ms1_peaks <- as.matrix(ms1_peaks[order(-ms1_peaks[, 'into']), ]) # arrange them by decreasing intensity; useful and main key in our strategy
    ms1_peaks.df <- as.data.frame(ms1_peaks)  } 

  if( "data.frame" %in% class(ms1_peaks) ){
    ms1_peaks.df <- ms1_peaks %>%
      arrange(desc(into))  }
  
  ms1_peaks.df$is_ion <- FALSE
  
  return(ms1_peaks.df)
}
