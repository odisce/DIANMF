#' Detect MS1 peaks using XCMS.
#'
#' @param rawData.onDiskMSnExp `OnDiskMSnExp` object for `onDisk` mode.
#' @param ... Additional parameters passed to \code{\link[xcms]{CentWaveParam}}.
#'
#' @return MS1 peaks `matrix`.
#' 
#' @export
#' 
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
#' 
#' @export
#' 
#' @examples
#' peaks_matrix <- matrix(c(150.5, 148.8, 155.2, 4.5, 4.9, 4.2, 3000,
#'                          200.3, 198.1, 205.7, 5.1, 5.5, 4.8, 4500,
#'                          180.7, 177.6, 185.0, 6.0, 6.4, 5.7, 5000,
#'                          210.9, 208.3, 215.5, 7.2, 7.6, 6.8, 3200),
#'                          nrow = 4, ncol = 7,
#'                          byrow = TRUE, dimnames = list(NULL, 
#'                          c("mz", "mzmin", "mzmax", "rt", "rtmax", "rtmin",
#'                           "into")))
#' 
#' prepare_ms1_peaks(ms1_peaks = peaks_matrix)
#' 
#' @importFrom dplyr arrange desc
#' 
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
