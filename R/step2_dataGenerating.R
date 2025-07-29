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
#' @param nscans number of minimal consecutive value to check
#'
#' @return `Logical` `TRUE` if the vector contains at least n
#'         consecutive non-zero values, `FALSE` otherwise.
has_n_consecutive_non_zero <- function(row, nscans = 4) {
  # Find consecutive non-zero values
  non_zero_streaks <- rle(row != 0)
  any(non_zero_streaks$lengths[non_zero_streaks$values] >= nscans)
}

#' Get scan time and align between ms events
#'
#' @inheritParams extract_xcms_peaks
#' @param rt_range retention time range.
#' @param sample_idx `NULL` to process all files else, `numeric` index of the sample.
#'
#' @return `list` of `data.table` `data.frame` objects: raw_dt and time_dic.
#'
#' @importFrom xcms filterRt filterFile
#' @import magrittr data.table
#' @import dplyr
#' @export
align_scans <- function(msexp, rt_range, sample_idx) {
  raw_dt <- xcms::filterFile(msexp, sample_idx) %>%
    suppressMessages()
  if (!is.null(rt_range)) {
    raw_dt <- xcms::filterRt(raw_dt, rt_range) %>%
      suppressMessages()
  }
  nev <- get_ms1_rtdiff(raw_dt) * 1.5
  raw_dt <- raw_dt %>%
    xcms::spectra(.) %>%
    Spectra::spectraData(.) %>%
    as.data.table()
  ## Normalize time
  time_dic <- raw_dt[, .(rtime, msLevel, isolationWindowTargetMz)] %>% unique()
  time_dic_ms1 <- time_dic[order(rtime), ][msLevel == 1, ]
  time_dic_ms1[, scan_norm := seq_len(.N)]
  ms2_rtime <- time_dic_ms1[
    ,
    {
      rt_ref <- rtime
      time_dic[rtime %between% (rt_ref + c(-nev, +nev)), .SD[which.min(abs(rtime - rt_ref)),], by = .(isolationWindowTargetMz)]
    },
    by = .(scan_norm)
  ]
  time_dic <- ms2_rtime[scan_norm %between% (range(scan_norm) + c(+1, -1)), ]
  time_dic[order(rtime), scan_norm := seq_len(.N), by = .(msLevel, isolationWindowTargetMz)]
  ## remove trailing scans
  scan_sel <- time_dic[, (range(scan_norm) + c(+1, -1)),]
  time_dic <- time_dic[scan_norm %between% scan_sel, ]
  time_dic[, scan_norm := scan_norm - 1]
  time_dic[]
}


# get_rawD_ntime <- function(
#   msexp,
#   rt_range,
#   sample_idx,
#   mslevel = 1
# ) {
#   raw_dt <- xcms::filterFile(msexp, sample_idx) %>%
#     xcms::filterRt(., rt_range)
#   if (!is.null(mslevel)) {
#     if (!mslevel %in% unique(msLevel(raw_dt))) {
#       stop(sprintf("mslevel %i not found in this file", mslevel))
#     } else {
#       raw_dt <- xcms::filterMsLevel(raw_dt, mslevel)
#     }
#   }
#   raw_dt <- get_spectra_values(raw_dt)
#   ## Normalize time
#   time_dic <- align_scans(
#     msexp = msexp,
#     rt_range = rt_range,
#     sample_idx = sample_idx
#   )
#   raw_dt <- merge(
#     raw_dt,
#     time_dic,
#     by = c("rtime", "msLevel", "isolationWindowTargetMz"),
#     allow.cartesian = TRUE
#   )
#   return(
#     list(
#       "raw_dt" = raw_dt,
#       "time_dic" = time_dic) # time dic no longer necessary because present in raw_dt...
#     )
# }


#' Build MS1 xics from peak list
#'
#' @param peaks_dt `data.table` `data.frame` targeted peaks in one rt_range.
#' @param rawdt `data.frame` from `DIANMF::get_rawD_ntime`.
#' @param method `string` to select a method for returning data:
#'   - `all`: return all points in range.
#'   - `max`: return max point by retention time.
#' @return `data.table` `data.frame` MS1 EICs.
#' @export
#'
#' @import data.table
build_XICs <- function(
  peaks_dt,
  rawdt,
  method = c("all", "max")[2]
) {
  grpcol <- intersect(c("peakid", "peakfull", "xic_label"), names(peaks_dt))

  
  if ("IsoWin" %in% names(peaks_dt)) {
    grid_dt <- peaks_dt[, .(msLevel, isolationWindowTargetMz = IsoWin)] %>% unique()
    key_join <- c("msLevel", "isolationWindowTargetMz")
  } else {
    grid_dt <- peaks_dt[, .(msLevel)] %>% unique()
    key_join <- c("msLevel")
  }

  output <- lapply(seq_len(grid_dt[, .N]), function(rowi) {
    tempi <- grid_dt[rowi,]
    x <- rawdt[tempi, on = key_join][,  .(start = mz, end = mz, mz, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy)]
    y <- peaks_dt[tempi, on = key_join][, .(start = mzmin, end = mzmax, xic_label)]
    
    if( anyNA(x$start) | anyNA(x$end) | anyNA(y$start) | anyNA(y$end)  ){
      
      return(NULL)
      
    }else {
      setkey(x, start, end)
      setkey(y, start, end)
      res <- foverlaps(x, y, type = "any", nomatch = NULL)
      if (method == "all") {
        out_i <- res[, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy, xic_label)]
      } else if (method == "max") {
        out_i <- res[,
                     .(
                       mz = median(mz, na.rm = TRUE),
                       intensity = max(intensity, na.rm = TRUE)
                     ),
                     by = .(rtime, msLevel, isolationWindowTargetMz, collisionEnergy, xic_label)
        ]
      } else {
        stop(sprintf("method argument not recognizes: %s", method))
      }
      return(out_i)
    }
  }) %>%
    rbindlist()
  return(output)
}

#' Build MS2 xics
#'
#' @inheritParams extract_xcms_peaks
#' @inheritParams get_spectra_values
#' @param ms2isowinL `logical` is `TRUE` to build the MS2 xics from raw data, else from xcms peaks. 
#' @param ppm `numeric` used for returning mzmin and mzmax values.
#'
#' @return `data.table` `data.frame` MS2 EICs.
#' @export
#'
#' @importFrom xcms filterFile filterRt filterMsLevel spectra isolationWindowTargetMz  
#' @importFrom Spectra filterIsolationWindow combineSpectra asDataFrame
#' @importFrom data.table as.data.table rbindlist
generate_peaklist <- function(
  msexp,
  mslevel = 1,
  ms2isowinL = TRUE,
  combineSpectra_arg,
  ppm = 6
) {
  if (ms2isowinL) {
    if (! 2 %in% mslevel) {
      stop("Can't extract isolation windows if mslevel argument = 1")
    }
    isowin <- get_isowin(msexp)
    if (length(isowin) <= 0) {
      stop("No isolation window found in the file")
    }
    export_grid <- data.table("ISOWIN" = isowin, MSLEVEL = 2)
  } else {
    export_grid <- data.table("ISOWIN" = as.numeric(NA), MSLEVEL = mslevel)
    export_grid[mslevel == 2, ISOWIN := "ALL"]
  }

  peaklist <- data.table()
  for (i in export_grid[, seq_len(.N)]) {
    isowini <- export_grid[i, ISOWIN]
    outi <- get_spectra_values(
      msexp = msexp,
      mslevel = export_grid[i, MSLEVEL],
      isowin = if (is.na(isowini) || isowini == "ALL") {NULL} else {isowini},
      combineSpectra_arg = combineSpectra_arg
    ) %>%
      {.[, .(IsoWin = isowini, msLevel, mz, rtime, intensity, isolationWindowTargetMz)]}
    peaklist <- rbind(peaklist, outi)
  }
  peaklist[, xic_label := paste0(seq_len(.N), "-", IsoWin, "-", msLevel)]
  peaklist[, c("mzmin", "mzmax") := {
    PpmRange(mz, ppm) %>% as.list()
  }, by = .(xic_label)]
  return(peaklist[])
}


# #' Build MS2 xics
# #'
# #' @inheritParams extract_xcms_peaks
# #' @inheritParams build_ms1XICS
# #' @param time_dic `data.frame` from `DIANMF::get_rawD_ntime`.
# #' @inheritParams get_rawD_ntime
# #'
# #' @return `data.table` `data.frame` MS2 EICs.
# #' @export
# #'
# #' @importFrom xcms filterFile filterRt filterMsLevel spectra isolationWindowTargetMz  
# #' @importFrom Spectra filterIsolationWindow combineSpectra asDataFrame
# #' @importFrom data.table as.data.table rbindlist
# build_XICs <- function(
#   msexp,
#   peaksdt
# ) {
#   ## Extract raw data

#   ## Extract XICs
#   xic_dt_ms2 <- MS2_peaklist[, {
#     mzrange <- PpmRange(mz, 7)
#     raw_dt[msLevel == 2 & mz %between% mzrange, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz,
#                                                   collisionEnergy)]
#   }, by = .(xic_label)]
#   xic_dt_ms2 <- xic_dt_ms2[xic_label %in% xic_dt_ms2[, .N, by = xic_label][N >= 4, xic_label]]
#   xic_dt_ms2 <- merge(
#     xic_dt_ms2,
#     time_dic,
#     by = c("rtime", "msLevel", "isolationWindowTargetMz"),
#     allow.cartesian = T
#   )
#   xic_dt_ms2[, xic_label := paste0(xic_label, "-", isolationWindowTargetMz)]
  
#   return(xic_dt_ms2)
# }

#' Create the mixed matrix
#' 
#' @param xicdt `data.table` of the xics build with `DIANMF::build_XICs()`
#' @inheritParams has_n_consecutive_non_zero
#'
#' @return `list` of the mixed matrices (ions x scans) with two levels:by
#'            - \[1\] mixed matrix with more than nscans consectuvies scans
#'            - \[2\] mixed matrix with less than nscans consectuvies scans
#' @export
build_mixed_matrix <- function(xicdt, nscans = 4) {
  mixeddt <- dcast(
    xicdt[order(scan_norm),],
    xic_label ~ scan_norm,
    value.var = "intensity",
    fun.aggregate = max,
    fill = 0
  )
  mixedmat <- as.matrix(mixeddt, rownames = 1)
  row_filter_ms1 <- apply(mixedmat, 1, has_n_consecutive_non_zero, n = nscans)
  
  # # put the eics of mslevel 1 at the upper part# it introduce some errors in the results!! I should check why
  # mixedmat <- mixedmat[grepl("-1$", rownames(mixedmat)), , drop = FALSE] %>%
  #   rbind(mixedmat[grepl("-2$", rownames(mixedmat)), , drop = FALSE])
  return(
    list(
      'mixedmat' = mixedmat[row_filter_ms1,],
      'mixedmat_deleted' = mixedmat[!row_filter_ms1,]
    )
  )
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
  row_filter_ms1 <- apply(ms1_mixedmat, 1, has_n_consecutive_non_zero)
  ms1_mixedmat <- ms1_mixedmat[row_filter_ms1, , drop = FALSE]
  
  ms1_mixedmat_deleted <- ms1_mixeddt[!row_filter_ms1, , drop = FALSE]
  
  ms1_infos <- xic_dt_ms1[msLevel == 1, ][, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]
  ms1_infos$isolationWindowLowerMz <- NA
  ms1_infos$isolationWindowUpperMz <- NA
  
  return(
    list(
      'ms1_mixedmat' = ms1_mixedmat,
      'ms1_mixedmat_deleted' = ms1_mixedmat_deleted,
      'ms1_infos' = ms1_infos)
    )
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
  row_filter_ms2 <- apply(ms2_mixedmat, 1, has_n_consecutive_non_zero)
  ms2_mixedmat <- ms2_mixedmat[row_filter_ms2, , drop = FALSE]
  ms2_infos <- xic_dt_ms2[
    ,
    .(mz = median(mz)),
    by = .(
      xic_label,
      msLevel,
      isolationWindowTargetMz
    )
  ]
  
  return(list(
    'ms2_mixedmat' = ms2_mixedmat,
    'ms2_infos' = ms2_infos
  ))
}

