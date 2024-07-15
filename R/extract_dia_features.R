# this script is written to process the raw files using xcms: step 1

#' get raw dia features from a mzML file
#'
#' @param files list_of_strings mzML files paths
#' @param d.out path results path
#' @param ppm.pd integer max mz devaition in consecutive scans
#' @param sn integer singl to noise ration
#' @param peakwidth integer range of the chromatogram peak width in seconds
#' @param mzdiff integer
#' @param nSlaves integer
#'
#' @importFrom xcms xcmsSet groups groupval
#' @importFrom BiocParallel SnowParam
#'
#' @return dia.feature extracted using xcms
#'
dia_features.f <- function(files, d.out, ppm.pd = 6, sn = 0, peakwidth = c(5,30), mzdiff = 0.0003, nSlaves){
  
  is.xcms3 <- FALSE
  if (packageVersion("xcms") >= '1.50.1'){
    is.xcms3 <- TRUE
    param <- SnowParam(workers = nSlaves, type = 'SOCK')
  }
  
  cat('Detect peaks with xcms ...\n')
  
  if (is.xcms3) {
    xset <- xcmsSet(files, method = 'centWave', ppm = ppm.pd, snthr = sn, peakwidth = peakwidth, mzdiff = mzdiff, BPPARAM = param)
  } else {
    xset <- xcmsSet(files,  method = 'centWave', ppm = ppm.pd, snthr = sn, peakwidth = peakwidth, mzdiff = mzdiff, nSlaves = nSlaves)
  }
  
  # save the xset object
  save(xset, file = paste0(d.out, '/xset-centWave.Rda'))
  # output the detected peak number
  sclassv <- as.character(sampclass(xset))
  for (i in 1:length( xset@filepaths)) {
    cat( i, xset@filepaths[ i], "\t [",  sclassv[ i], "]", sep = "")
    cat( "  --> ", length( which( xset@peaks[, "sample"] == i)), " Features. \n")
  }
  
  return(xset)
}

#-------------------------------------------------------------------------------

#' get raw dia features from a mzML file
#'
#' @param d.in mzML files paths
#' @param d.out path results path
#' @param nSlaves integer
#' @param peakwidth integer range of the chromatogram peak width in seconds
#' @param mzdiff integer
#' @param sn integer singl to noise ration
#' @param ppm.pd integer max mz devaition in consecutive scans
#' @param rerun boolean
#' @param correlation boolean
#' 
#' @importFrom xcms xcmsSet groups groupval
#' @importFrom BiocParallel SnowParam
#'
#' @return extract all dia features of all mzML files
#'
extract_dia_features.f <- function(d.in,
                                      d.out, nSlaves = 6,
                                      # peak detection setup
                                      peakwidth = c(5, 30), mzdiff = 0.001, sn = 0, ppm.pd = 6,
                                      rerun = FALSE,
                                      # to process multiple samples
                                      correlation = FALSE ) {
  wd0 <- getwd();
  files <- list.files(d.in, recursive = TRUE, full.names = TRUE, pattern = '(?i)mzML$')
  
  if(correlation == TRUE ){
    d.out <- paste0(wd0, '/', d.out);
    if (!dir.exists(d.out)) {
      dir.create(d.out, recursive = TRUE)
    };
    
    nSlaves <- min(detectCores() - 1, nSlaves, length(file))
    bpparam <- SnowParam(workers = nSlaves, type = 'SOCK')
    CatLineSeparator('Detecting and aligning features ...')
    fn.skip <- paste0(d.out ,'/dia.feature.RData')
    
    if ((!rerun) & file.exists(fn.skip)) {
      cat('using existing results:', fn.skip, '...\n')
      load(fn.skip)
    } else {
      dia.feature <- dia_features.f(files, d.out, ppm.pd = 6, sn = 0, peakwidth = c(5,30), mzdiff = 0.0003, nSlaves)
    }
  } else {
    for(i in 1:length(files)){
      one_file <- files[[i]]
      d.out <- paste0(wd0, '/', d.out, '/', file_path_sans_ext(basename(one_file)))
      if (!dir.exists(d.out)) {
        dir.create(d.out, recursive = TRUE)
      };
      
      nSlaves <- min(detectCores() - 1, nSlaves, length(file))
      bpparam <- SnowParam(workers = nSlaves, type = 'SOCK')
      CatLineSeparator('Detecting and aligning features ...')
      fn.skip <- paste0( d.out ,'/dia.feature.RData')
      
      if ((!rerun) & file.exists(fn.skip)) {
        cat('using existing results:', fn.skip, '...\n')
        load(fn.skip)
      } else {
        dia.feature <- dia_features.f(files, d.out, ppm.pd = 6, sn = 0, peakwidth = c(5,30), mzdiff = 0.0003, nSlaves)
      }
    }
  }
  
  return(dia.feature)
}

#-------------------------------------------------------------------------------
