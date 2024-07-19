
#' Detect EICs peaks
#'
#' @param rawData.onDiskMSnExp 
#' @param d.out path results path
#' @param ppm integer max mz devaition in consecutive scans
#' @param peakwidth integer range of the chromatogram peak width in seconds
#' @param snthresh integer singl to noise ration
#' @param prefilter integer
#' @param mzCenterFun string method
#' @param integrate integer
#' @param mzdiff integer
#'
#' @importFrom xcms findChromPeaks
#'
#' @return rawData.XCMSnExp
#'
detect_EICs.f <- function(rawData.onDiskMSnExp, d.out = NULL,
                      ppm = 6, peakwidth = c(6,60), snthresh = 1,
                      prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                      integrate = 2, mzdiff = -0.001, noise = 0, firstBaselineCheck = FALSE
                      ){
  
  cat('Detect peaks with xcms ...\n')
  # create centwave parameter object
  cwp <- xcms::CentWaveParam(ppm = ppm,
                             peakwidth = peakwidth,
                             snthresh = snthresh,
                             prefilter = prefilter,
                             mzCenterFun = mzCenterFun,
                             integrate = integrate,
                             mzdiff = mzdiff,
                             noise = noise,
                             firstBaselineCheck = firstBaselineCheck,
                             );
  # detect EIC peaks
  rawData.XCMSnExp <- findChromPeaks(rawData.onDiskMSnExp, param = cwp);
    
  # save the object
  if (!is.null(d.out)) {
    if (!dir.exists(d.out)) {
      dir.create(d.out)
    }
    save(rawData.XCMSnExp, file = paste0(d.out, '/rawData.XCMSnExp'))
  }
  # print the detected peak number
  # ...............................
  
  return(rawData.XCMSnExp)
}


# extract_rawData.f <- function(d.in, d.out, nSlaves = 6,
#                               # peak detection setup
#                               ppm = 6, snthresh = 1, peakwidth = c(6, 60), mzdiff = -0.001,
#                               rerun = FALSE,
#                               # to process multiple samples
#                               correlation = FALSE ){
#   
#   wd0 <- getwd(); # the working directory path 
#   files <- list.files(d.in, recursive = TRUE, full.names = TRUE, pattern = '(?i)mzML$'); # list of all mzML files in the d.in directory
#   
#   if(correlation == TRUE ){  # to process all files at once
#     d.out <- paste0(wd0, '/', d.out);  # create the output directory to save results
#     if (!dir.exists(d.out)) {  
#       dir.create(d.out, recursive = TRUE)
#     };
#     
#     rawData.onDiskMSnExp <- readMSData(one_file, mode = "onDisk");  # read mzML as OnDiskMSnExp data from MSnbase package
#     
#     fn.skip <- paste0(d.out ,'/rawData.XCMSnExp'); # check if the raw_MSnbase_data already exist; if not create them
#     if ((!rerun) & file.exists(fn.skip)) {
#       cat('using existing results:', fn.skip, '...\n')
#       load(fn.skip)
#     } else {
#       nSlaves <- min(detectCores() - 1, nSlaves, length(file)) 
#       rawData.XCMSnExp <- detect_EICs.f(rawData.onDiskMSnExp, d.out, ppm = ppm, snthresh = snthresh, peakwidth = peakwidth, mzdiff = mzdiff, nSlaves)
#     }
#   } else {  # to process every file separately
#     for(i in 1:length(files)){
#       one_file <- files[[i]]
#       d.out <- paste0(wd0, '/', d.out, '/', file_path_sans_ext(basename(one_file)))
#       if (!dir.exists(d.out)) {
#         dir.create(d.out, recursive = TRUE)
#       };
#       
#       rawData.onDiskMSnExp <- readMSData(one_file, mode = "onDisk");  # read mzML as OnDiskMSnExp data from MSnbase package
#       
#       fn.skip <- paste0( d.out ,'/rawData.XCMSnExp')  # check if the EIC peaks already detected
#       if ((!rerun) & file.exists(fn.skip)) {
#         cat('using existing results:', fn.skip, '...\n')
#         load(fn.skip)
#       } else {
#         nSlaves <- min(detectCores() - 1, nSlaves, length(file))
#         rawData.XCMSnExp <- detect_EICs.f(rawData.onDiskMSnExp, d.out, ppm = ppm, snthresh = snthresh, peakwidth = peakwidth, mzdiff = mzdiff, nSlaves)
#       }
#     }
#   }
#   
#   return(raw_MSnbase_data)
# }
