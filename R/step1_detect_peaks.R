#' Detect MS1 peaks using XCMS.
#'
#' @param sequence_table A data.frame with at least [mzml_path, class, InjectionOrder].
#' @param params A list to set the parameters from xcms pipeline steps.
#'               List names should be xcms parameters methods (ex: `list("CentWaveParam" = CentWaveParam()`)
#' @return MsExperiment objects
#' 
#' @export
#' 
#' @import MsExperiment MsFeatures xcms data.table magrittr
detect_xcms_peaks <- function(
  sequence_table,
  params = list(
    "CentWaveParam" = xcms::CentWaveParam(),
    "MergeNeighboringPeaksParam" = xcms::MergeNeighboringPeaksParam(),
    "ObiwarpParam" = xcms::ObiwarpParam(),
    "PeakDensityParam" = xcms::PeakDensityParam(sampleGroups = NA),
    "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
  )
){
  sequence_table <- data.table::as.data.table(sequence_table)
  sequence_table <- sequence_table[order(InjectionOrder), ]
  xcms_obj <- MsExperiment::readMsExperiment(spectraFiles = sequence_table$mzml_path, sampleData = sequence_table)
  if (length(params$PeakDensityParam@sampleGroups) == 1 && is.na(params$PeakDensityParam@sampleGroups)) {
    params$PeakDensityParam@sampleGroups <- sampleData(xcms_obj)$class
  }
  xcms_obj_peaks <- xcms::findChromPeaks(
    xcms_obj,
    param = params$CentWaveParam
  ) %>% 
    refineChromPeaks(., param = params$MergeNeighboringPeaksParam, msLevel = 1L)
  if (length(xcms::fileNames(xcms_obj_peaks)) > 1) {
    xcms_obj_peaks <- adjustRtime(xcms_obj_peaks, param = params$ObiwarpParam)
  }
  output <- groupChromPeaks(xcms_obj_peaks, param = params$PeakDensityParam)
  
  if (length(xcms::fileNames(xcms_obj_peaks)) > 1) {
    output <- fillChromPeaks(output, param = params$ChromPeakAreaParam)
  }
  return(output)
}

#' Detect MS1 peaks using XCMS.
#'
#' @param MsExperiment.obj `MsExperiment` object obtained from xcms or with `DIANMF::detect_xcms_peaks()`.
#' @param sample_nb Index of the sample to extract peaks from.
#'
#' @return MS1 peaks `matrix`.
#' 
#' @export
#' 
#' @importFrom xcms findChromPeaks CentWaveParam chromPeaks
#' @import magrittr MsExperiment
extract_xcms_peaks <- function(MsExperiment.obj, sample_nb = 1) {
  filterFile(MsExperiment.obj, sample_nb) %>%
    xcms::chromPeaks(.) %>%
    return()
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
  
  ms1_peaks.df$is_ion <- 0
  
  return(ms1_peaks.df)
}
