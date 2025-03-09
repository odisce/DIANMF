# Step2: extract MS1 and MS2 info and matrices

#' Create a mz range
#'
#' @param ref `numeric` mz of the peak.
#' @param ppm.n `numeric` defining the m/z tolerated deviation in parts per million (ppm).
#'
#' @return `numeric` mz range.
PpmRange <-  function(ref, ppm.n) {
  dev <- ppm.n * 1e-6
  ref + (c(-1, 1) * ref * dev)
}


#' Check if the vector has at least 4 consecutive non-zero values
#'
#' @param row vector of `numeric` to evaluate.
#'
#' @return `Logical` `TRUE` if the vector contains at least 4  consecutive non-zero values, `FALSE` otherwise.
has_four_consecutive_non_zero <- function(row) {
  # Find consecutive non-zero values
  non_zero_streaks <- rle(row != 0)
  any(non_zero_streaks$lengths[non_zero_streaks$values] >= 4)
}


#' Extract MS1 raw data
#'
#' @inheritParams extract_xcms_peaks
#' @param rt_range retention time range.
#' @param sample_idx `numeric` sample index.
#'
#' @return `list` of `data.table` `data.frame` objects: raw_dt and time_dic.
#' @export
#' 
#' @importFrom xcms filterRt filterFile
#' @import magrittr data.table
#' @import dplyr
get_rawD_ntime <- function(msexp, rt_range, sample_idx){
  
  raw_dt <- xcms::filterRt(msexp, rt_range) %>%
    xcms::filterFile(., sample_idx) %>%
    get_spectra_values()
  
  ## Normalize time
  time_dic <- raw_dt[, .(rtime, msLevel, isolationWindowTargetMz)] %>% unique()
  time_dic_ms1 <- time_dic[order(rtime), ][msLevel == 1, ]
  time_dic_ms1[, scan_norm := seq_len(.N)]
  ms2_rtime <- time_dic_ms1[, {
    rt_ref <- rtime
    time_dic[order(abs(rtime - rt_ref)), ][, lapply(.SD, head, 1), by = .(isolationWindowTargetMz)]
  }, by = .(scan_norm)]
  time_dic <- ms2_rtime[scan_norm %between% (range(scan_norm) + c(+1, -1)), ]
  time_dic[order(rtime), scan_norm := seq_len(.N), by = .(msLevel, isolationWindowTargetMz)]
  
  raw_dt <- merge(
    raw_dt,
    time_dic,
    by = c("rtime", "msLevel", "isolationWindowTargetMz"),
    allow.cartesian = T  
    )
  
  return(list(
    'raw_dt' = raw_dt,
    'time_dic' = time_dic  ))
}


#' Build MS1 xics from peak list
#'
#' @param peaks_i `data.table` `data.frame` targeted peaks in one rt_range.
#' @param raw_dt `data.frame` from `DIANMF::get_rawD_ntime`.
#'
#' @return `data.table` `data.frame` MS1 EICs.
#' @export
#' 
#' @import data.table 
#' @importFrom dplyr between
build_ms1XICS <- function(peaks_i, raw_dt){
  
  xic_dt_ms1 <- peaks_i[, {
    mzrange <- c(mzmin, mzmax)
    raw_dt[msLevel == 1 & mz %between% mzrange, .(mz, scan_norm, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy)]
  }, by = .(peakid, peakfull)]
  xic_dt_ms1[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
  xic_dt_ms1 <- xic_dt_ms1[xic_label %in% xic_dt_ms1[, .N, by = xic_label][N >= 4, xic_label]]  # filter MS1 eics
  
  return(xic_dt_ms1)
}

#' Build MS2 xics
#'
#' @inheritParams extract_xcms_peaks
#' @inheritParams build_ms1XICS
#' @param time_dic `data.frame` from `DIANMF::get_rawD_ntime`.
#' @inheritParams get_rawD_ntime
#' @param MS2_ISOEACHL `logical` is TRUE to build the MS2 xics from raw data, else from xcms peaks. 
#'
#' @return `data.table` `data.frame` MS2 EICs.
#' @export
#'
#' @importFrom xcms filterFile filterRt filterMsLevel spectra isolationWindowTargetMz  
#' @importFrom Spectra filterIsolationWindow combineSpectra asDataFrame
#' @importFrom data.table as.data.table rbindlist
build_ms2XICs <- function(msexp, raw_dt, time_dic, rt_range, MS2_ISOEACHL = T){
  
  # build MS2 xics from raw data
  if (MS2_ISOEACHL) {
    isowin <- msexp %>%
      xcms::filterFile(., 1) %>%
      xcms::filterRt(., rt_range) %>%
      xcms::filterMsLevel(., 2L) %>%
      xcms::spectra() %>%
      xcms::isolationWindowTargetMz() %>%
      as.numeric() %>% unique()
    MS2_peaklist_l <- list()
    for (n in seq_len(length(isowin))) {
      i <- isowin[[n]]
      MS2_peaklist_l[[n]] <- xcms::filterRt(msexp, rt_range + c(+1, -1)) %>%
        xcms::filterFile(., 1) %>%
        xcms::filterMsLevel(., 2) %>%
        xcms::spectra() %>%
        Spectra::filterIsolationWindow(i) %>%
        Spectra::combineSpectra(ppm = 3, tolerance = 0.005, minProp = 0.001) %>%
        Spectra::asDataFrame() %>%
        data.table::as.data.table() %>%
        {.[, .(msLevel, mz, rtime, intensity, isolationWindowTargetMz)]}
    }
    MS2_peaklist <- data.table::rbindlist(MS2_peaklist_l, idcol = "IsoWin")
    MS2_peaklist[, xic_label := paste0(seq_len(.N), "-", IsoWin, "-", 2)]
  } else {
    MS2_peaklist <- xcms::filterRt(msexp, rt_range) %>%
      xcms::filterFile(., 1) %>%
      xcms::filterMsLevel(., 2) %>%
      xcms::spectra() %>%
      Spectra::combineSpectra(ppm = 5, tolerance = 0.008, minProp = 0.001) %>%  # (ppm = 5, tolerance = 0.008, minProp = 0.01)
      Spectra::asDataFrame() %>%
      data.table::as.data.table() %>%
      {.[, .(msLevel, mz, rtime, intensity)]}
    MS2_peaklist[, xic_label := paste0(seq_len(.N), "-", 2)]
  }
  
  # Extract MS2 XICs
  xic_dt_ms2 <- MS2_peaklist[, {
    mzrange <- PpmRange(mz, 7)
    raw_dt[msLevel == 2 & mz %between% mzrange, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz,
                                                  collisionEnergy)]
  }, by = .(xic_label)]
  xic_dt_ms2 <- xic_dt_ms2[xic_label %in% xic_dt_ms2[, .N, by = xic_label][N >= 4, xic_label]]
  xic_dt_ms2 <- merge(
    xic_dt_ms2,
    time_dic,
    by = c("rtime", "msLevel", "isolationWindowTargetMz"),
    allow.cartesian = T
  )
  xic_dt_ms2[, xic_label := paste0(xic_label, "-", isolationWindowTargetMz)]
  
  return(xic_dt_ms2)
}


#' Generate MS1 data and info
#' 
#' @param xic_dt_ms1 `data.table` `data.frame` object obtained from xcms or with `DIANMF::build_ms1XICs()`.
#'
#' @return `list` of extracted MS1 data (MS1 matrix, deleted matrix, MS1 info).
#' @export
ms1Info <- function(xic_dt_ms1){
  ms1_mixeddt <- dcast(xic_dt_ms1[msLevel == 1, ], xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
  row_names <- ms1_mixeddt[, 1]
  ms1_mixedmat <- ms1_mixeddt <- as.matrix(ms1_mixeddt[, -1])
  rownames(ms1_mixedmat) <- rownames(ms1_mixeddt) <- t(row_names)
  row_filter_ms1 <- apply(ms1_mixedmat, 1, has_four_consecutive_non_zero)
  ms1_mixedmat <- ms1_mixedmat[row_filter_ms1, , drop = FALSE]
  
  ms1_mixedmat_deleted <- ms1_mixeddt[!row_filter_ms1, , drop = FALSE]
  
  ms1_infos <- xic_dt_ms1[msLevel == 1, ][, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]
  ms1_infos$isolationWindowLowerMz <- NA
  ms1_infos$isolationWindowUpperMz <- NA
  
  return(list(
    'ms1_mixedmat' = ms1_mixedmat,
    'ms1_mixedmat_deleted' = ms1_mixedmat_deleted,
    'ms1_infos' = ms1_infos))
}


#' Generate MS2 data and info
#'
#' @param xic_dt_ms2  `data.table` `data.frame` object obtained from xcms or with `DIANMF::build_ms2XICs()`.
#'
#' @return `list` of extracted MS2 data (MS2 matrix, MS2 info).
#' @export 
ms2Info <- function(xic_dt_ms2){
  ms2_mixeddt <- dcast(xic_dt_ms2, xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
  row_names <- ms2_mixeddt[, 1]
  ms2_mixedmat <- as.matrix(ms2_mixeddt[, -1])
  rownames(ms2_mixedmat) <- t(row_names)
  # ms2_mixedmat <- as.matrix(ms2_mixeddt, rownames = TRUE)
  row_filter_ms2 <- apply(ms2_mixedmat, 1, has_four_consecutive_non_zero)
  ms2_mixedmat <- ms2_mixedmat[row_filter_ms2, , drop = FALSE]
  ms2_infos <- xic_dt_ms2[, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz)]
  
  return(list(
    'ms2_mixedmat' = ms2_mixedmat,
    'ms2_infos' = ms2_infos
  ))
}

