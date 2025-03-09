#' This is data to be included in my package
#'
#' @name data_example
#' @docType data
#' @author Diana Karaki \email{diana.karaki@cea.fr}
#' @keywords data
NULL

# @references \url{data_blah.com}

#' Get mzml test paths
#'
#' @return Vector of mzml paths
#'
#' @export
get_test_mzml <- function() {
  mzml_dir <- system.file("/inst/extdata", package = "DIANMF")
  mzml_files <- list.files(mzml_dir, pattern = "\\.mzml$", full.names = TRUE)
  return(mzml_files)
}

#' Load test sequence
#'
#' @return `data.table` with mzML paths
#' @import data.table
#'
#' @export
get_test_sequence <- function() {
  data.table(
    "mzml_path" = get_test_mzml(),
    "class" = c("A", "B"),
    "InjectionOrder" = c(1, 2)
  )
}

#' Perform peak picking on test data
#'
#' @return `XCMSnExp` object containing detected peaks from both mzML files
#' @import xcms data.table MsExperiment
#'
#' @export
get_test_peaks <- function(){
  sequence_table <- get_test_sequence()
  
  sequence_table %>%
    detect_xcms_peaks(
      .,
      params = list(
        "CentWaveParam" = xcms::CentWaveParam(ppm = 10, peakwidth = c(3, 10), snthresh = 2, firstBaselineCheck = FALSE),
        "MergeNeighboringPeaksParam" = xcms::MergeNeighboringPeaksParam(),
        "ObiwarpParam" = xcms::ObiwarpParam(),
        "PeakDensityParam" = xcms::PeakDensityParam(sampleGroups = sequence_table$class),
        "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
      )
    )
}

