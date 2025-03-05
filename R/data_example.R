#' This is data to be included in my package
#'
#' @name data_example
#' @docType data
#' @author Diana Karaki \email{diana.karaki@cea.fr}
#' @keywords data
NULL

# @references \url{data_blah.com}

#' Get mzml test path
#'
#' @return mzml path
#'
#' @export
get_test_mzml <- function() {
  mzml_path <- system.file("/inst/extdata/test_data.mzml", package = "DIANMF")
}

#' Load test data
#'
#' @return MsExperiment objects
#' @import MsExperiment
#'
#' @export
get_test_data <- function() {
  MsExperiment::readMsExperiment(get_test_mzml())
}

#' Load test sequence
#'
#' @return MsExperiment objects
#' @import data.table
#'
#' @export
get_test_sequence <- function() {
  data.table(
    "mzml_path" = get_test_mzml(),
    "class" = "A",
    "InjectionOrder" = 1
  )
}

#' Perform peakpicking on test data
#'
#' @return MsExperiment objects
#' @import data.table xcms
#'
#' @export
get_test_peaks <- function() {
  get_test_sequence() %>%
    detect_xcms_peaks(
      .,
      params = list(
        "CentWaveParam" = xcms::CentWaveParam(ppm = 10, peakwidth = c(3,10), snthresh = 2, firstBaselineCheck = FALSE),
        "MergeNeighboringPeaksParam" = xcms::MergeNeighboringPeaksParam(),
        "ObiwarpParam" = xcms::ObiwarpParam(),
        "PeakDensityParam" = xcms::PeakDensityParam(sampleGroups = NA),
        "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
      )
  )
}
