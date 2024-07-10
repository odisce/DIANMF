# this script let us choose a peak that is related to a precursor and extract the mixed matrix of the ms1 and ms2 data that I want

#--------------------------------------------- Let's choose a peak

#' extract rt, scans, eics of a psecific MS1 peak
#'
#' @param pg.raw list raw peaks
#' @param fp.smp file sample
#' @param idx.pg index of the peak
#' @param scanTimes_mzRange MS1 & MS2 scans information (scans times and spectrum for every scan)
#' @param rt.info rt information for every MS1 peak
#' @param fp.swath dictionary learning file
#' @param int.filter minimum intensity to filter
#' @param ppm.ms2.mtp ppm value to pick ms2 peaks
#' @param dia.feature dia features
#'
#' @import DecoMetDIA
#'
#' @return rt, scans indices, apex info of ms1 & ms2 peaks and eics raw, extended, and smooth
#'
info_all <- function(pg.raw, fp.smp, idx.pg, scanTimes_mzRange, rt.info, fp.swath, int.filter, ppm.ms2.mtp, dia.feature){

  idx.smp <- which(dia.feature@filepaths %in% fp.smp)

  idx.ms1.raw <- idx.ms1 <- which(scanTimes_mzRange$scantime.ms1 >= rt.info[[1]][idx.pg, 'rtmin'] &
                                    scanTimes_mzRange$scantime.ms1 <= rt.info[[1]][idx.pg, 'rtmax']) # get the indices of the ms1 scans where this peak may belong

  idx.ms1.ext.raw <- idx.ms1.ext <- which(scanTimes_mzRange$scantime.ms1 >= rt.info[[1]][idx.pg, 'rtmin.ext'] &
                                            scanTimes_mzRange$scantime.ms1 <= rt.info[[1]][idx.pg, 'rtmax.ext']) # get the indices of the ms1 scans where the extended peak may belong

  mz <- pg.raw[[idx.pg]][idx.smp, 'mz'] # mz of this peak
  mz.max <- PpmRange(ref = mz, ppm = 25)[2] # create mz range and take mzmax
  mz.range <- as.numeric(pg.raw[[idx.pg]][idx.smp, c('mzmin', 'mzmax')]) # extract the true mz range of this peak

  info.swath <- GetSWATHinfo(fp = fp.swath) # load the ms2 windows (nSWATH = nb of ms2 windows)
  scanidx.match <- MatchScanidx(fp.smp, info.swath) # match the ms2 scans to there ms1 scans
  idx.swath <- GetIntervalIndex(mz, info.swath$window) # get the index of swath window where this ms1 peak of this mz value belongs
  idx.ms2 <- unname(sapply(scanidx.match[as.character(idx.ms1)], `[`, idx.swath)) # extract the ms2 scans(related to this peak) that belongs to idx.swath
  idx.ms2.ext <- unname(sapply(scanidx.match[as.character(idx.ms1.ext)], `[`, idx.swath)) # extract the ms2 scans(related to the extended peak) that belongs to idx.swath

  eic.ms1.pre <- extractEIC(scanTimes_mzRange$scan.ms1[idx.ms1.ext],
                            matrix(mz.range, ncol = 2),
                            mz)[[1]] # list of length 2(mz column, intensity column); ms1 EIC of this peak
  peaks.ms1.ext <- cbind('scanidx' = idx.ms1.ext,
                         'rt' = scanTimes_mzRange$scantime.ms1[idx.ms1.ext],
                         'mz' = eic.ms1.pre$mz,
                         'intensity' = eic.ms1.pre$intensity) # transform eic.ms1.pre as a matrix

  is.keep <- idx.ms1.ext %in% idx.ms1
  peak.ms1 <- peaks.ms1.ext[is.keep, , drop = FALSE] # keep just the ms1 points that are related to this peak

  peak.ms1.ext.smooth.raw <- peak.ms1.ext.smooth <- SmoothLoess(peaks.ms1.ext, span = 0.1) # smooth the ms1 points
  peak.ms1.smooth.raw <- peak.ms1.smooth <- peak.ms1.ext.smooth[is.keep, , drop = FALSE] # keep just the smoothed ms1 points that are related to this peak

  idx.na <- which(is.na(idx.ms2.ext)) # get the indices of the null scans
  if (length(idx.na) > 0) {
    idx.ms2.ext <- idx.ms2.ext[-idx.na]
    if (length(idx.ms2.ext) == 0) {
      return(NULL)
    }
    idx.ms1.na <- idx.ms1[idx.ms1 == idx.ms1.ext[idx.na]]
    if (length(idx.ms1.na) > 0) {
      peak.ms1.ext.smooth <- peak.ms1.ext.smooth[-idx.ms1.na, , drop = FALSE]
      peaks.ms1.ext <- peaks.ms1.ext[-idx.ms1.na, , drop = FALSE]
      idx.ms1 <- idx.ms1.raw[-idx.ms1.na]
      idx.ms1.ext <-idx.ms1.ext[-idx.ms1.na]
    }
  } # remove the null scans

  if (all(peak.ms1.smooth[, 'intensity.s'] == 0)) {
    return(NULL)
  } # if there is no peaks ...

  # index of the apex scan without extension
  idx.apex <- idx.ms1[which.max(peak.ms1.smooth[,'intensity.s'])]

  idx.apex.ms1 <- which(idx.ms1.ext == idx.ms1[which.max(peak.ms1.smooth[,'intensity.s'])]) # get the index of the ms1 scan where the apex of this ms1 peak belongs

  rt.apex.ms1 <- peak.ms1.ext.smooth[idx.apex.ms1, 'rt'] # get the rt of the apex of this ms1 peak
  if (nrow(peak.ms1.smooth) < 4) {
    return(NULL)
  } # if the ms1 peak contains less than 4 points we will not consider it as a peak

  spec.exp.ext <- scanTimes_mzRange$scan.ms2[idx.ms2.ext] # get the spectra of the ms2 scans
  spec.apex <- spec.exp.ext[[idx.apex.ms1]] # get the ms2 spectrum where the ms1 apex belongs
  idx.ms2.intfilter <- which(spec.apex[,2] >= int.filter &
                               spec.apex[, 1] <= mz.max) # filter the ions in spec.apex
  mz.ms2 <- spec.exp.ext[[idx.apex.ms1]][idx.ms2.intfilter, 1] # mz values
  if (length(mz.ms2) == 0) {
    return(NULL)
  }
  mz.ms2.range <- t(sapply(mz.ms2, function(mz) {
    PpmRange(mz, ppm.ms2.mtp)
  }))

  eic.ms2.pre <- extractEIC(spec.exp.ext, mz.ms2.range, mz.ms2)
  ms2.eic.ext <- lapply(eic.ms2.pre, function(spec) {
    cbind('scanidx' = idx.ms2.ext,
          'rt' = scanTimes_mzRange$scantime.ms2[idx.ms2.ext],
          do.call(cbind, spec))
  })

  range.pk.ms1 <- range(which(is.keep))
  info.pk.ms1 <- data.frame('lb' = range.pk.ms1[1],
                            'apex' = idx.apex.ms1,
                            'apex.t' = rt.apex.ms1,
                            'rb' = range.pk.ms1[2],
                            'bl' = pg.raw[[idx.pg]][idx.smp, ]$maxo/pg.raw[[idx.pg]][idx.smp, ]$sn,
                            'simp' = FALSE)
  rownames(info.pk.ms1) <- 0

  return(list(
    "idx.ms2.ext" = idx.ms2.ext,
    "eic.ms1.pre" = eic.ms1.pre,
    "peaks.ms1.ext" = peaks.ms1.ext,
    "peak.ms1" = peak.ms1,
    "peak.ms1.ext.smooth"  = peak.ms1.ext.smooth,
    "peak.ms1.smooth" = peak.ms1.smooth,
    "idx.apex.ms1" = idx.apex.ms1,
    "rt.apex.ms1" = rt.apex.ms1,
    "eic.ms2.pre" = eic.ms2.pre,
    "ms2.eic.ext" = ms2.eic.ext,
    "info.pk.ms1" = info.pk.ms1,
    "spec.apex" = spec.apex,
    "idx.ms2" = idx.ms2,
    "idx.ms1.raw" = idx.ms1.raw,
    'idx.apex' = idx.apex ))
}

#-------------------------------------------------------------------------------

#' extract the mixed MS1 data matrix
#'
#' @param pg.raw list raw peaks
#' @param fp.smp string file path
#' @param idx.pg numeric index of the MS1 peak
#' @param scanTimes_mzRange list MS1 & MS2 scans, rt of the scans and the mz range
#' @param rt.info list retention time information raw and extended for MS1 peak
#' @param fp.swath string SWATH file path
#' @param int.filter numeric minimum intensity accepted
#' @param ppm.ms2.mtp numeric used to extract the mz range
#' @param dia.feature list contains peaks information of the file
#'
#' @import DecoMetDIA
#'
#' @return list of MS1 information related the MS1 peak idx.pg; the MS1 matrix, rt and mz values of this peak
#'
get_part1_matrix_X <- function(pg.raw, fp.smp, idx.pg, scanTimes_mzRange, rt.info, fp.swath, int.filter, ppm.ms2.mtp, dia.feature){

  idx.smp <- which(dia.feature@filepaths %in% fp.smp)

  idx.ms1.raw <- idx.ms1 <- which(scanTimes_mzRange$scantime.ms1 >= rt.info[[1]][idx.pg, 'rtmin'] &
                                    scanTimes_mzRange$scantime.ms1 <= rt.info[[1]][idx.pg, 'rtmax']) # get the indices of the ms1 scans where this peak may belong

  mz <- pg.raw[[idx.pg]][idx.smp, 'mz'] # mz of this peak
  mz.range <- as.numeric(pg.raw[[idx.pg]][idx.smp, c('mzmin', 'mzmax')]) # extract the true mz range of this peak

  eic.ms1.pre <- extractEIC(scanTimes_mzRange$scan.ms1[idx.ms1],
                            matrix(mz.range, ncol = 2),
                            mz)[[1]] # list of length 2(mz column, intensity column); ms1 EIC of this peak
  peaks.ms1 <- cbind('scanidx' = idx.ms1,
                     'rt' = scanTimes_mzRange$scantime.ms1[idx.ms1],
                     'mz' = eic.ms1.pre$mz,
                     'intensity' = eic.ms1.pre$intensity) # transform eic.ms1.pre as a matrix

  idx.apex.ms1 <- which(idx.ms1 == idx.ms1[which.max(peaks.ms1[,'intensity'])]) # get the index of the ms1 scan where the apex of this ms1 peak belongs

  info.swath <- GetSWATHinfo(fp = fp.swath) # load the ms2 windows (nSWATH = nb of ms2 windows)
  idx.swath <- GetIntervalIndex(mz, info.swath$window) # get the index of swath window where this ms1 peak of this mz value belongs
  isol_window_mz_range <- as.numeric(info.swath$window[idx.swath, ]) # get the mz range of this isolation window (I will use this info to filter the ions in the ms1 apex spectra )

  # extract the ms1 information ---------------------------------------------------------------------------------------------

  # get the ms1 scans related to the ms1 peak
  spec.exp_ms1 <- scanTimes_mzRange$scan.ms1[idx.ms1]

  # get the ms1 scan that contains the ms1 peak apex
  spec.apex_ms1 <- spec.exp_ms1[[idx.apex.ms1]]

  # filter this spectrum
  idx.ms1.filter <- which(spec.apex_ms1[ ,2] >= int.filter &
                            spec.apex_ms1[ ,1] >= isol_window_mz_range[1] &
                            spec.apex_ms1[ ,1] <= isol_window_mz_range[2] ) # filter the ions in spec.apex_ms1 by intensity and isolation window mz range

  # get the mz_values remains in the spectrum after filtering
  mz.ms1 <- spec.exp_ms1[[idx.apex.ms1]][idx.ms1.filter, 1] # mz values

  # create their mz range
  mz.ms1.range <- t(sapply(mz.ms1, function(mz) {
    PpmRange(mz, ppm.ms2.mtp)
  }))

  # get the eics of the previous mz values from the ms1 scans related to the ms1 peak
  eic.ms1.pre_ms1 <- extractEIC(spec.exp_ms1, mz.ms1.range, mz.ms1)

  ms1.eic_ms1 <- lapply(eic.ms1.pre_ms1, function(spec) {
    cbind('scanidx' = idx.ms1,
          'rt' = scanTimes_mzRange$scantime.ms1[idx.ms1],
          do.call(cbind, spec))
  })

  rt.apex.ms1 <- peaks.ms1[idx.apex.ms1, 'rt'] # get the rt of the apex of this ms1 peak
  if (nrow(peaks.ms1) < 4) {
    return(NULL)
  } # if the ms1 peak contains less than 4 points we will not consider it as a peak

  # Done ---------------------------------------------------------------------------------------------

  # creating the mixed matrix the part related to ms1

  # extract the intensity vectors from eic.ms1.pre_ms1
  ms1_intensities <- lapply(1:length(eic.ms1.pre_ms1), function(i){
    inten <- eic.ms1.pre_ms1[[i]]$intensity
  })
  ms1_intensity.m <- do.call(rbind, ms1_intensities)
  rownames(ms1_intensity.m) <- paste0('MS1_',round(mz.ms1,4))

  # extract the intensity vector from the ms1 EIC of the ms1 peak
  # eic_ms1_intensity <- eic.ms1.pre$intensity
  #
  # X <- rbind(eic_ms1_intensity, ms1_intensity.m)
  # rownames(X)[1] <- 'MS1'
  X <- ms1_intensity.m
  #---------------------------------------------- Done

  # extract the retention time vector (it is common btw all the ms1 eics related to the apex)
  ms1_rt <- ms1.eic_ms1[[1]][, 'rt']

  # extract the scans idx
  scans_ids <- ms1.eic_ms1[[1]][ ,'scanidx']

  # extract the mz vector
  ms1_mz <- lapply(1:length(ms1.eic_ms1), function(i){
    mz_value <- mean(ms1.eic_ms1[[i]][, 'mz'])

  })

  # ms1_mz <- c(mz, ms1_mz)  # at the begining of the list I added new level that contains the ms1 peak mz  or bzed mz of the apex???!!! answer it is the same
  ms1_mz <- ms1_mz
  return(list( 'X.m' = X,
               'rt.v' = ms1_rt,
               'mz.l' = ms1_mz))
}

#-------------------------------------------------------------------------------

#' extract the mixed MS2 data matrix, while trying to decrease the shifting between the MS1 and MS2 data
#'
#' @param pg.raw list raw peaks
#' @param fp.smp string file path
#' @param idx.pg numeric index of the MS1 peak
#' @param scanTimes_mzRange list MS1 & MS2 scans, rt of the scans and the mz range
#' @param rt.info list retention time information raw and extended for MS1 peak
#' @param fp.swath string SWATH file path
#' @param int.filter numeric minimum intensity accepted
#' @param ppm.ms2.mtp numeric used to extract the mz range
#' @param dia.feature list contains peaks information of the file
#'
#' notes: when the peak may belon to scan index 1 I cant extend, so no retention time alignment
#' also some MS1 peaks wasnt fragmented in the MS2 steps so I cant extract there MS2 data; they doesnt exists
#'
#' @import DecoMetDIA
#'
#' @return list of MS2 information related the MS1 peak idx.pg; the MS2 matrix, rt and mz values of this peak
#'
get_part2_matrix_X <- function(pg.raw, fp.smp, idx.pg, scanTimes_mzRange, rt.info, fp.swath, int.filter, ppm.ms2.mtp, dia.feature){

  idx.smp <- which(dia.feature@filepaths %in% fp.smp)  # sample index

  idx.ms1 <- which(scanTimes_mzRange$scantime.ms1 >= rt.info[[1]][idx.pg, 'rtmin'] &
                     scanTimes_mzRange$scantime.ms1 <= rt.info[[1]][idx.pg, 'rtmax']) # get the indices of the ms1 scans where this peak may belong

  idx.ms1.semi.ext <- c(idx.ms1[1] -1, idx.ms1, idx.ms1[length(idx.ms1)] + 1) # partially extend idx.ms1 (1 scan from left & 1 scan from right)

  # idx.ms1.semi.ext <- idx.ms1.semi.ext[idx.ms1.semi.ext != 0] # exclude the scan 0 since it doesn't exist !!!!!!!!!!!!!!!!!!!!!!!!!

  mz <- pg.raw[[idx.pg]][idx.smp, 'mz'] # mz of this peak
  mz.range <- as.numeric(pg.raw[[idx.pg]][idx.smp, c('mzmin', 'mzmax')]) # extract the true mz range of this peak

  info.swath <- GetSWATHinfo(fp = fp.swath) # load the ms2 windows (nSWATH = nb of ms2 windows)

  # if the peak wasn't fragmented in the MS2 step, exist; there is no MS2 data related to it !!!!!!!!!!!!!!!!!!!!!!!!!
  if(mz < info.swath$window[1,1] | mz > info.swath$window[nrow(info.swath$window),2]) { return(NULL) };

  idx.swath <- GetIntervalIndex(mz, info.swath$window) # get the index of swath window where this ms1 peak of this mz value belongs

  scanidx.match <- MatchScanidx(fp.smp, info.swath) # match the ms2 scans to there ms1 scans
  idx.ms2.ext <- unname(sapply(scanidx.match[as.character(idx.ms1.semi.ext)], `[`, idx.swath)) # the ms2 scans(related to this extended peak) from idx.swath

  idx.ms2 <- unname(sapply(scanidx.match[as.character(idx.ms1)], `[`, idx.swath)) # the ms2 scans related to this MS1 peak from idx.swath


  eic.ms1.pre <- extractEIC(scanTimes_mzRange$scan.ms1[idx.ms1.semi.ext],
                            matrix(mz.range, ncol = 2),
                            mz)[[1]] # list of length 2(mz column, intensity column); ms1 extended EIC of this peak

  peaks.ms1.ext <- cbind('scanidx' = idx.ms1.semi.ext,
                         'rt' = scanTimes_mzRange$scantime.ms1[idx.ms1.semi.ext],
                         'mz' = eic.ms1.pre$mz,
                         'intensity' = eic.ms1.pre$intensity) # transform eic.ms1.pre as a matrix and add some info

  is.keep <- idx.ms1.semi.ext %in% idx.ms1
  peak.ms1 <- peaks.ms1.ext[is.keep, , drop = FALSE] # keep just the ms1 points that are related to this peak not for the extended one

  idx.na <- which(is.na(idx.ms2.ext)) # get the indices of the null scans
  if (length(idx.na) > 0) {
    idx.ms2.ext <- idx.ms2.ext[-idx.na]
    if (length(idx.ms2.ext) == 0) {
      return(NULL)
    }
  } # remove the null scans

  if (all(peak.ms1[, 'intensity'] == 0)) {
    return(NULL)
  } # if there is no peaks ...

  if (nrow(peak.ms1) < 4) {
    return(NULL)
  } # if the ms1 peak contains less than 4 points we will not consider it as a peak


  if(idx.swath <= 5 | 1 %in% idx.ms1){  # !!!!!!!!!!!!!!!!!!!!!!!!! also I should add the condition if the scan index is 10 in the 10th isolation window
    # this is not exactly the good condition, itshould be if isolation window is 1 and idx.ms1 ==1 & isolation window is 10 and idx.ms1 == 10

    idx.apex.ms1 <- which(idx.ms1 == idx.ms1[which.max(peak.ms1[,'intensity'])]) # index of the ms1 scan that contains the apex

    spec.exp <- scanTimes_mzRange$scan.ms2[idx.ms2] # the spectra of the ms2 scans
    spec.apex <- spec.exp[[idx.apex.ms1]] # the ms2 spectrum where the ms1 apex belongs
    idx.ms2.intfilter <- which(spec.apex[,2] >= int.filter) # filter the ions in spec.apex
    mz.ms2 <- spec.exp[[idx.apex.ms1]][idx.ms2.intfilter, 1] # mz values that we want to keep from the ms2 apex spectrum
    if (length(mz.ms2) == 0) {
      return(NULL)
    }  # if the ms2 apex spectrum is empty, return NULL

    mz.ms2.range <- t(sapply(mz.ms2, function(mz) {
      PpmRange(mz, ppm.ms2.mtp)
    }))  # make a mz range for every me in the apex MS2 spectrum

    eic.ms2.pre <- extractEIC(spec.exp, mz.ms2.range, mz.ms2)  # MS2 eics related to this MS1 peak
    ms2.eic.ext <- lapply(eic.ms2.pre, function(spec) {
      cbind('scanidx' = idx.ms2,
            'rt' = scanTimes_mzRange$scantime.ms2[idx.ms2],
            do.call(cbind, spec))
    })  # transform every eic in eic.ms2.pre for a matrix and add some info

  } else {

    idx.apex.ms1 <- which(idx.ms1.semi.ext == idx.ms1[which.max(peak.ms1[,'intensity'])]) # get the index of the ms1 scan where the apex of this ms1 peak belongs

    spec.exp <- scanTimes_mzRange$scan.ms2[idx.ms2.ext] # get the spectra of the ms2 scans
    spec.apex <- spec.exp[[idx.apex.ms1]] # get the ms2 spectrum where the ms1 apex belongs
    idx.ms2.intfilter <- which(spec.apex[,2] >= int.filter) # filter the ions in spec.apex
    mz.ms2 <- spec.exp[[idx.apex.ms1]][idx.ms2.intfilter, 1] # mz values
    if (length(mz.ms2) == 0) {
      return(NULL)
    }
    mz.ms2.range <- t(sapply(mz.ms2, function(mz) {
      PpmRange(mz, ppm.ms2.mtp)
    }))

    eic.ms2.pre <- extractEIC(spec.exp, mz.ms2.range, mz.ms2)

    ms2.eic.ext <- lapply(eic.ms2.pre, function(spec) {
      cbind('scanidx' = idx.ms2.ext,
            'rt' = scanTimes_mzRange$scantime.ms2[idx.ms2.ext],
            do.call(cbind, spec))
    })

    eic.ms2.pre <- lapply(1:length(eic.ms2.pre), function(i){
      res <- lapply(eic.ms2.pre[[i]], function(x) head(x, -2))
    })  # MS2 eics related to this MS1 peak

    ms2.eic.ext <- lapply(1:length(ms2.eic.ext), function(i){
      res <- ms2.eic.ext[[i]][1:(nrow(ms2.eic.ext[[i]]) - 2), ]
    }) # transform every eic in eic.ms2.pre for a matrix and add some info

  }

  # creating the MS2 mixed matrix
  # 1- non extended part
  ms2_intensities <- lapply(1:length(eic.ms2.pre), function(i){
    inten <- eic.ms2.pre[[i]]$intensity
  }) # extract the intensity vectors from eic.ms2.pre

  ms2_intensity.m <- do.call(rbind, ms2_intensities) # put them in a matrix
  rownames(ms2_intensity.m) <- paste0('MS2_',round(mz.ms2,4)) # assign every row (=a intensity vector of an eic) to its eic name
  X <- ms2_intensity.m

  #---------------------------------------------- Done

  # extract the retention time vector (it is common between all the ms2 eics related to the apex)
  ms2_rt <- ms2.eic.ext[[1]][, 'rt']

  # extract the scans idx
  scans_ids <- ms2.eic.ext[[1]][ ,'scanidx']

  # extract the mz vector
  ms2_mz <- lapply(1:length(ms2.eic.ext), function(i){
    mz_value <- mean(ms2.eic.ext[[i]][, 'mz'])
  })

  return(list( 'X.m' = X,
               'rt.v' = ms2_rt,
               'mz.l' = ms2_mz))

}

#-------------------------------------------------------------------------------
