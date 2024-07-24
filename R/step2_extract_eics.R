# extract the isolation windows of SWATH data

#' Title
#'
#' @param rawData.onDiskMSnExp 
#'
#' @return isolation windows
#' @export
#' @import MSnbase
#' 
isolationWindows.range <- function(rawData.onDiskMSnExp){
  
  ms2_spectra <- MSnbase::filterMsLevel(rawData.onDiskMSnExp, msLevel = 2) # filter MS2 scans
  # Extract isolation window lower and upper m/z values
  isolation_windows <- data.frame(
    lowerMz = unique(MSnbase::isolationWindowLowerMz(ms2_spectra)),
    upperMz = unique(MSnbase::isolationWindowUpperMz(ms2_spectra))
  )
  
  return(isolation_windows)
}

#' create a mz range
#'
#' @param ref 
#' @param ppm 
#'
#' @return mz range of a specific mz value
#' @export
#'
PpmRange <-  function(ref, ppm) {
  dev <- ppm * 1e-6
  ref + (c(-1, 1) * ref * dev)
}

# to be fixed: I want to filter the ions of the MS1 apex spectra, by deleteing the ions of intensity lower than ms1_int_filter

#' Title
#'
#' @param spectra_list 
#' @param ppm 
#' @param apex_index 
#' @param rt_index 
#' @param mz_range 
#'
#' @return extracted eics
#' @export
#' @import data.table
#' @import magrittr
#'
extract_eics <- function(spectra_list, ppm = 7, apex_index, rt_index = TRUE, mz_range = TRUE) {
  
  full_table <- data.table::rbindlist(spectra_list, idcol = "spectra_index")
  if (!is.null(mz_range)) {
    full_table <- full_table[mz %between% mz_range,]
  }
  
  target_table <- full_table[spectra_index == apex_index,]
  target_table[, index := seq_len(.N)]
  target_table <- target_table[, {
    out <- PpmRange(mz, ppm)
    .(start = min(out), end = max(out))
  }, by = .(index)]
  full_table[, start := mz][, end := mz]
  setkey(target_table, start, end)
  setkey(full_table, start, end)
  
  output <- data.table::foverlaps(
    target_table,
    full_table
  )
  output <- merge(output, unique(full_table[, .(spectra_index)]), by = "spectra_index", all = TRUE)
  
  output[, meanmz := mean(mz) %>% sprintf("%.4f", .), by = .(index)]
  output <- output[order(-intensity), lapply(.SD, head, 1), by = .(spectra_index, meanmz)]
  
  # require(ggplot2)
  # ggplot(output, aes(spectra_index, intensity, color = index, group = index)) +
  #   geom_line() +
  #   theme_bw()
  mixed_table <- dcast(output, meanmz ~ spectra_index, value.var = "intensity")
  # max_ind <- apply(mixed_table, 1, max, na.rm = T) %>% which.max()
  # mixed_table <- dcast(output, meanmz ~ spectra_index, value.var = "mz")
  # mixed_table[max_ind,]
  mixed_matrix <- as.matrix(mixed_table, rownames = 1)
  if (!is.null(rt_index)) {
    colnames(mixed_matrix) <- rt_index[as.integer(colnames(mixed_matrix))] %>% sprintf("%.3f", .)
  }
  mixed_matrix[is.na(mixed_matrix)] <- 0
  return(mixed_matrix)
}
# spectra_list <- spec.exp_ms1
# ppm <- 7
# apex_index <- idx.apex.ms1
# mz_range <- isol_window_mz_range
# mz_range <- NULL
# rt_index <- spec.exp_rt