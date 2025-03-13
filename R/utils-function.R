#' Extract the isolation windows of SWATH data.
#'
#' @param x MsExperiment object
#' @return `data.frame` SWATH isolation windows.
#' 
#' @export
#' 
#' @import magrittr xcms data.table
fun_get_rawdata <- function(x) {
  spctr <- spectra(x)
  mz_values <- lapply(spctr, mz)
  int_values <- lapply(spctr, intensity)
  rt_values <- lapply(spctr, rtime)
  file_value <- lapply(spctr, fromFile)

  output <- data.table(
    "file" = rep(unlist(file_value), lengths(mz_values)),
    "scan" = rep(1:length(mz_values), lengths(mz_values)),
    "rt" = rep(unlist(rt_values), lengths(mz_values)),
    "mz" = unlist(mz_values),
    "i" = unlist(int_values)
  )
  setkey(output, file, scan, rt, mz)
  return(output)
}

#' Get unique msLevel from MsExperiment object
#'
#' @param msexp MsExperiment object
#' @return an integer vector with available MSLevels
#' 
#' @export
#' 
#' @import magrittr MsExperiment
get_mslevels <- function(msexp) {
  Spectra::uniqueMsLevels(MsExperiment::spectra(msexp)) %>%
    return()
}

#' Get spectra value in MsExperiment
#'
#' @param msexp MsExperiment object
#' @param mslevel (optional) MSLevel to subset
#' @param isowin (optional) isolation window to subset
#' @param combineSpectra_arg (optional) combineSpectra argument as a list
#' @return a data.table containing spectra infos 
#' @import MsExperiment
#' @importFrom data.table as.data.table
#' @importFrom Spectra asDataFrame filterIsolationWindow combineSpectra
#' @export
get_spectra_values <- function(msexp, mslevel = NULL, isowin = NULL, combineSpectra_arg = NULL) {
  spec_obj <- MsExperiment::spectra(msexp)
  if (!is.null(mslevel)) {
    spec_obj <- Spectra::filterMsLevel(spec_obj, mslevel)
  }
  if (!is.null(isowin)) {
    spec_obj <- Spectra::filterIsolationWindow(spec_obj, isowin)
  }
  if (!is.null(combineSpectra_arg)) {
    spec_obj <- do.call(Spectra::combineSpectra, args = c(list("x" = spec_obj), combineSpectra_arg)) %>%
      suppressMessages()
  }
  spec_obj %>%
    Spectra::asDataFrame() %>%
    data.table::as.data.table() %>%
    return()
}

#' Get spectra index in MsExperiment
#'
#' @param msexp MsExperiment object
#' @return a data.table with scan infos
#' 
#' @export
#' 
#' @import MsExperiment Spectra data.table
get_spectra_index <- function(msexp) {
  Spectra::spectraData(spectra(msexp)) %>%
    data.table::as.data.table() %>%
    return()
}

#' Get rt diff between 2 MS1 scans
#'
#' @param msexp MsExperiment object
#' @return a data.table containing spectra infos
#' @importFrom xcms filterMsLevel rtime
#' @export
get_ms1_rtdiff <- function(msexp) {
  msexp %>%
    xcms::filterMsLevel(., 1L) %>%
    xcms::rtime() %>%
    {.[1:2]} %>%
    diff() %>%
    abs()
}

#' Get isolation windows
#'
#' @param msexp MsExperiment object
#' @return a numeric vector
#' @importFrom xcms spectra isolationWindowTargetMz
#' @export
get_isowin <- function(msexp) {
  msexp %>%
    xcms::spectra() %>%
    xcms::filterMsLevel(., 2) %>%
    xcms::isolationWindowTargetMz() %>%
    as.numeric() %>%
    unique() %>%
    na.omit() %>%
    return()
}
