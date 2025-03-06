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
#' @return `MsExperiment` objects
#' @import MsExperiment
#'
#' @export
get_test_data <- function() {
  MsExperiment::readMsExperiment(get_test_mzml())
}

#' Load test sequence
#'
#' @return `MsExperiment` objects
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

#' Perform peak picking on test data
#'
#' @return `MsExperiment` objects
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


# params_ls <- list(
#   "CentWaveParam" = CentWaveParam(
#     ppm = 6,
#     peakwidth = c(6, 30),
#     snthresh = 0,
#     prefilter = c(5, 4000),
#     mzCenterFun = "wMeanApex3",
#     integrate = 2,
#     mzdiff = -0.0003,
#     noise = 2000,
#     firstBaselineCheck = FALSE
#   ),
#   "MergeNeighboringPeaksParam" = MergeNeighboringPeaksParam(
#     expandRt = 2,
#     expandMz = 0,
#     ppm = 1,
#     minProp = 0.75
#   ),
#   "ObiwarpParam" = ObiwarpParam(
#     binSize = 0.05
#   ),
#   "PeakDensityParam" = PeakDensityParam(
#     sampleGroups = NA,
#     bw = 15,
#     minFraction = 0,
#     minSamples = 1,
#     binSize = 0.005,
#     ppm = 5,
#     maxFeatures = 500
#   ),
#   "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
# )
