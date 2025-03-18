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
  spec_obj <- setBackend(spec_obj, MsBackendMemory())
  if (!is.null(mslevel)) {
    spec_obj <- Spectra::filterMsLevel(spec_obj, mslevel) %>% suppressMessages()
  }
  if (!is.null(isowin)) {
    spec_obj <- Spectra::filterIsolationWindow(spec_obj, isowin) %>% suppressMessages()
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
    suppressMessages() %>%
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
    suppressMessages() %>%
    xcms::isolationWindowTargetMz() %>%
    as.numeric() %>%
    unique() %>%
    na.omit() %>%
    return()
}

#' Get isolation windows
#'
#' @param feature_dt feature summary
#' @param dt database to match (mandatory rt and mz columns)
#' @param rttol rt tolerance (in seconds)
#' @param ppm ppm tolerance for m/Z search
#' @return a data.table
#' @import data.table magrittr
#' @export
search_features <- function(
  feature_dt,
  dt,
  rttol = 5,
  ppm = 5
) {
    dtin <- copy(dt)
    dtin[, TEMPID := seq_len(.N)]
    dtin[, c("mzmin", "mzmax", "rtmin", "rtmax") := {
      tpi <- PpmRange(mz, ppm)
      tprt <- (rt + c(-rttol, +rttol))
      .(tpi[1], tpi[2], tprt[1], tprt[2])
    }, by = .(TEMPID)]
    setnames(dtin, c("mz", "rt"), c("mz_db", "rt_db"))
    feat_in <- copy(feature_dt)
    feat_in <- feat_in[, .(mz = median(mz), rt = median(rt)), by = .(featureid)]
    resmz <- data.table()
    if (!is.null(ppm)) {
      ## Filter mz
      x <- feat_in[,  .(start = mz, end = mz, featureid)] %>% unique()
      y <- dtin[, .(start = mzmin, end = mzmax, TEMPID)]
      setkey(x, start, end)
      setkey(y, start, end)
      resmz <- foverlaps(x, y, type = "any", nomatch = NULL)[, .(TEMPID, featureid, match = "mz")]
    }
    if (!is.null(rttol)) {
      ## Filter rt
      x <- feat_in[,  .(start = rt, end = rt, featureid)] %>% unique()
      y <- dtin[, .(start = rtmin, end = rtmax, TEMPID)]
      setkey(x, start, end)
      setkey(y, start, end)
      resrt <- foverlaps(x, y, type = "any", nomatch = NULL)[, .(TEMPID, featureid, match = "rt")]
    }
    out <- rbind(resmz, resrt)
    out_match <- out[, .(matchL = .N == 2), by = .(featureid, TEMPID)][matchL == TRUE,]
    out_match %>%
      merge(
        .,
        feat_in,
        by = "featureid"
      ) %>%
      merge(
        .,
        dtin,
        by = "TEMPID"
      ) %>%
      return()
}