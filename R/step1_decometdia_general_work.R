# this script use decometdia functions the do the general work; i.e. to get the ms1 & ms2 scans information and rt info for every ms1 peak

#' extract ms1 & ms2 scans information
#'
#' @param fp.smp file sample
#' @param dia.feature dia features
#'
#' @import xcms
#'
#' @return ms1 & ms2 scans information (scans times and spectrum for every scan)
scan_time.l<- function( fp.smp, dia.feature){

  idx.smp <- which(dia.feature@filepaths %in% fp.smp) # get index of the sample file for these dia.feautures
  if (length(idx.smp) == 0) {
    if (file.exists(fp.smp)) {
      idx.smp <- which(basename(dia.feature@filepaths) %in% basename(fp.smp))
    } else {
      stop(paste0("file not exist: \n  ", fp.smp))
    }
  }

  xr <- xcmsRaw(fp.smp, includeMSn = TRUE) # transform fp.smp data to xcmsRaw object so we can use them
  scantime.ms1 <- xr@scantime # extract the ms1 scans times
  mz.range <- xr@mzrange # extract the mz range
  scan.ms1 <- lapply(seq_along(xr@scantime) , function(idx) {
    getScan(xr, idx)
  })  # for every ms1 scan (ms1 mass scan), get the mz and intensity values; for every ms1 scan get the spectrum

  xr.ms2 <- msn2xcmsRaw(xr) # copy the MS2 data from xr to xcmsRaw object slots so we can use them
  scantime.ms2 <- xr.ms2@scantime # extract the ms2 scans times
  scan.ms2 <- lapply(seq_along(xr.ms2@scantime) , function(idx) {
    getScan(xr.ms2, idx)
  }) # for every ms2 scan (ms2 mass scan), get the mz and intensity values; for every ms2 scan get the spectrum

  return(list(
    "scantime.ms1" = scantime.ms1, # ms1 scans times
    "mz.range" = mz.range, # mz range
    "scan.ms1" = scan.ms1, # ms1 spectra
    "scantime.ms2" = scantime.ms2, # ms1 scans times
    "scan.ms2" = scan.ms2 # ms2 spectra
  ))
}



#' extract rt information
#'
#' @param pg.raw ms1 raw peaks
#' @param rt.extend rt value, to extend the peak rt range
#' @param scanTimes_mzRange  ms1 & ms2 scans information (scans times and spectrum for every scan)
#'
#' @return rt information for every ms1 peak
rt_info <- function(pg.raw, rt.extend, scanTimes_mzRange){

  seq.pg <- seq_along(pg.raw) # create seq 1:length(ms1 peaks)

  rt.info <- do.call(rbind, lapply(seq.pg, function(idx.pg) {
    pk <- pg.raw[[idx.pg]]
    data.frame('idx' = idx.pg, pk[pk$sample == 1, c('rt', 'rtmin', 'rtmax')])
  })) # data-frame, every raw contains info about 1 ms1 peak (its: idx, rt, rtmin, rtmax)

  if (!is.null(rt.extend)) {
    rt.info$rtmin.ext <- rt.info$rt - rt.extend * (rt.info$rtmax - rt.info$rtmin)
    rt.info$rtmax.ext <- rt.info$rt + rt.extend * (rt.info$rtmax - rt.info$rtmin)
    rt.info$rtmin.ext[rt.info$rtmin.ext < 0] <- 0
    rt.info$rtmax.ext[rt.info$rtmax.ext > max(scanTimes_mzRange$scantime.ms1)] <- max(scanTimes_mzRange$scantime.ms1)
  } else {
    rt.info$rtmin.ext <- rt.info$rtmin
    rt.info$rtmax.ext <- rt.info$rtmax
  } # data-frame, every raw contains info about 1 ms1 peak (its: idx, rt, rtmin, rtmax, rtmin.ext, rtmax.ext); the extended rt range for every ms1 peak was added

  return(list(
    "rt.info" = rt.info
  ))
}

#--------------------------------------------- General work finished

