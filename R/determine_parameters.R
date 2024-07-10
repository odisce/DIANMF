#' determine xcms parameters: for peak detection and alignment peaks from different samples
#'
#' parameters to detect peaks
#' @param ppm integer
#' @param sn integer
#' @param peakwidth integer
#' @param mzdiff integer
#'
#' @param nSlaves integer
#' @param rt.extend integer
#'
#' parameters to align peaks of multiple samples
#' @param bw integer
#' @param mzwid integer
#' @param minsamp integer
#' @param minfrac integer
#'
#' @return list of parameters for xcms
#' @export
#'
xcms_parameters <- function(ppm = 6, sn = 0, peakwidth = c(5, 30), mzdiff = 0.0003,
                            nSlaves = 6, rt.extend = 2.5,
                            bw = 5, mzwid = 0.25, minsamp = 0, minfrac = 0.5){
  return(list(
    'ppm' = ppm,
    'sn' = sn,
    'peakwidth' = peakwidth,
    'mzdiff' = mzdiff,
    'nSlaves' = nSlaves,
    'rt.extend' = rt.extend,
    'bw' = bw,
    'mzwid' = mzwid,
    'minsamp' = minsamp,
    'minfrac' = minfrac  ))
}

#-------------------------------------------------------------------------------

#' specify the parameters for MS1 & MS2 data extraction for every MS1 peak
#'
#' @param ppm.ms2.mtp integer
#' @param int.filter integer
#'
#' @return list of parameters to extract mixed matrices and there information
#' @export
#'
extract_ms1_ms2_data_parameters <- function(ppm.ms2.mtp = 15, int.filter = 50){
  return(list(
    'ppm.ms2.mtp' = ppm.ms2.mtp,
    'int.filter' = int.filter ))
}

#-------------------------------------------------------------------------------

#' specify the factorization algorithm parameters
#'
#' @param maximumIteration integer
#' @param maxFBIteration integer
#' @param toleranceFB integer
#' @param useTranspose Boolean
#'
#' @return list of parameters to run the nGMCAs algorithm
#' @export
#'
factorization_parameters <- function(maximumIteration = 50, maxFBIteration = 10, toleranceFB = 0.0001, useTranspose = TRUE){
  return(list(
    'maximumIteration' = maximumIteration,
    'maxFBIteration' = maxFBIteration,
    'toleranceFB' = toleranceFB,
    'useTranspose' = useTranspose ))
}

#-------------------------------------------------------------------------------
