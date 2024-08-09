#' Detect EICs peaks
#'
#' @param rawData.onDiskMSnExp raw MSnbase data of mzMl file
#' @param ppm integer max mz devaition in consecutive scans
#' @param peakwidth integer range of the chromatogram peak width in seconds
#' @param snthresh integer singl to noise ration
#' @param prefilter integer
#' @param mzCenterFun string method
#' @param integrate integer
#' @param mzdiff integer
#' @param noise numeric
#' @param firstBaselineCheck Boolean 
#'
#' @import MSnbase
#' @importFrom xcms findChromPeaks CentWaveParam
#'
#' @return rawData.XCMSnExp
#' @export
detect_peaks_by_xcms <- function(rawData.onDiskMSnExp,
                      ppm = 6, peakwidth = c(6,60), snthresh = 1,
                      prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                      integrate = 2, mzdiff = -0.001, noise = 0, firstBaselineCheck = FALSE
                      ){
  
  cat('Detect peaks with xcms ...\n')
  # create centwave parameter object
  cwp <- xcms::CentWaveParam(ppm = ppm,
                             peakwidth = peakwidth,
                             snthresh = snthresh,
                             prefilter = prefilter,
                             mzCenterFun = mzCenterFun,
                             integrate = integrate,
                             mzdiff = mzdiff,
                             noise = noise,
                             firstBaselineCheck = firstBaselineCheck,
                             );
  # detect EIC peaks
  rawData.XCMSnExp <- xcms::findChromPeaks(rawData.onDiskMSnExp, param = cwp);
    
  # # save the object
  # if (!is.null(d.out)) {
  #   if (!dir.exists(d.out)) {
  #     dir.create(d.out)
  #   }
  #   save(rawData.XCMSnExp, file = paste0(d.out, '/rawData.XCMSnExp'))
  # }
  
  return(rawData.XCMSnExp)
}


#' Transform the ms1_peaks R object to data.frame, and the rows are ordered in decreasing order of the 'into' column.
#' 
#' @param ms1_peaks matrix or data.frame contains all MS1 peaks 
#'
#' @return data.frame MS1 peaks
#' @export
#' @importFrom dplyr arrange desc
#' @import magrittr
prepare_ms1_peaks <- function(ms1_peaks){
  
  if( !("matrix" %in% class(ms1_peaks)) & !("data.frame" %in% class(ms1_peaks)) ){
    print("Error; the ms1_peaks should be a matrix or a data.frame R objects")
    return(NULL)
  };
  
  if( "matrix"  %in% class(ms1_peaks) ){
    ms1_peaks <- as.matrix(ms1_peaks[order(-ms1_peaks[, 'into']), ]) # arrange them by decreasing intensity; useful and main key in our strategy
    ms1_peaks.df <- as.data.frame(ms1_peaks)
  } 

  if( "data.frame" %in% class(ms1_peaks) ){
    ms1_peaks.df <- ms1_peaks %>%
      arrange(desc(into))
  }
  
  ms1_peaks.df$is_ion <- FALSE
  
  return(ms1_peaks.df)
}
