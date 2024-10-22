#' Extract the isolation windows of SWATH data.
#'
#' @inheritParams detect_peaks_by_xcms
#'
#' @return `data.frame` SWATH isolation windows.
#' 
#' @export
#' 
#' @importFrom MSnbase isolationWindowLowerMz isolationWindowUpperMz filterMsLevel
isolationWindows.range <- function(rawData.onDiskMSnExp){
  
  ms2_spectra <- MSnbase::filterMsLevel(rawData.onDiskMSnExp, msLevel = 2) # filter MS2 scans
  # Extract isolation window lower and upper m/z values
  isolation_windows <- data.frame(
    lowerMz = unique(MSnbase::isolationWindowLowerMz(ms2_spectra)),
    upperMz = unique(MSnbase::isolationWindowUpperMz(ms2_spectra))
  )
  
  return(isolation_windows)
}


#' Create a mz range.
#'
#' @param ref `numeric` mz of the peak.
#' @inheritParams extract_ms_matrix.f
#'
#' @return `numeric` mz range of a specific mz value.
PpmRange <-  function(ref, ppm.n) {
  dev <- ppm.n * 1e-6
  ref + (c(-1, 1) * ref * dev)
}


#' Extract the MS1 or MS2 eics from a list of spectra.
#'
#' @param spectra_list `list` of spectra.
#' @inheritParams extract_ms_matrix.f
#' @param apex_index `integer` spectrum index which contains the peak apex.
#' @param rt_index `Logical` if `TRUE` use the real retention time axis of the peak, else `FALSE`.
#' @param mz_range `Logical` if `TRUE` and extracting MS1 data , if `TRUE` and extracting MS2 data , else `FALSE`.
#'
#' @return `matrix` contains extracted eic in every row.
#' 
#' @importFrom data.table rbindlist foverlaps setkey merge.data.table dcast %between%
#' @import magrittr
#' @importFrom stats end start
#' @importFrom utils head
extract_eics <- function(spectra_list, ppm.n = 7, apex_index, rt_index = TRUE, mz_range = TRUE) {
  
  full_table <- data.table::rbindlist(spectra_list, idcol = "spectra_index")
  if (!is.null(mz_range)) {
    full_table <- full_table[mz %between% mz_range,]
  }
  
  target_table <- full_table[spectra_index == apex_index,]
  target_table[, index := seq_len(.N)]
  target_table <- target_table[, {
    out <- PpmRange(mz, ppm.n)
    .(start = min(out), end = max(out))
  }, by = .(index)]
  full_table[, start := mz][, end := mz]
  data.table::setkey(target_table, start, end)
  data.table::setkey(full_table, start, end)
  
  output <- data.table::foverlaps(
    target_table,
    full_table
  )
  output <- data.table::merge.data.table(output, unique(full_table[, .(spectra_index)]), by = "spectra_index", all = TRUE)
  
  output[, meanmz := mean(mz) %>% sprintf("%.4f", .), by = .(index)]
  output <- output[order(-intensity), lapply(.SD, head, 1), by = .(spectra_index, meanmz)]
  
  mixed_table <- data.table::dcast(output, meanmz ~ spectra_index, value.var = "intensity")
  mixed_matrix <- as.matrix(mixed_table, rownames = 1)
  if (!is.null(rt_index)) {
    colnames(mixed_matrix) <- rt_index[as.integer(colnames(mixed_matrix))] %>% sprintf("%.3f", .)
  }
  mixed_matrix[is.na(mixed_matrix)] <- 0
  return(mixed_matrix)
}