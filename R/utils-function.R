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
#' @return a data.table containing spectra infos
#' 
#' @export
#' 
#' @import MsExperiment Spectra data.table
get_spectra_values <- function(msexp) {
  Spectra::asDataFrame(MsExperiment::spectra(msexp)) %>%
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