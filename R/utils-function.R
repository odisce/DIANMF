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