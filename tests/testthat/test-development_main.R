skip("development script")

test_that("test main on barbier data", {
  require(xcms)
  require(MSnbase)
  require(ggplot2)
  require(dplyr)

  devtools::load_all()
  
  file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
  rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");
  fast_sample <- FALSE
  
  if (fast_sample) {
    sample_file <- MSnbase::filterRt(rawData.onDiskMSnExp, c(0,2*60))
  } else {
    sample_file <- rawData.onDiskMSnExp
  }
  eics_peaks <- detect_EICs.f(
    rawData.onDiskMSnExp = sample_file,
    d.out = NULL,
    ppm = 6,
    peakwidth = c(6,60), 
    snthresh = 1,
    prefilter = c(5,4000), 
    mzCenterFun = "wMeanApex3",
    integrate = 2,
    mzdiff = -0.001, 
    noise = 0, 
    firstBaselineCheck = FALSE
  )
  eics_peaks.mat <- xcms::chromPeaks(eics_peaks)   # extract the peaks detected by xcms
  eics_peaks.mat <- as.matrix(eics_peaks.mat[order(-eics_peaks.mat[, 'into']), ]) # arrange them by decreasing intensity
  
  nmf_parameters.l <-  list(
    'maximumIteration' = 10,
    'maxFBIteration' = 5,
    'toleranceFB' = 1e-5,
    'useTranspose' = TRUE,
    'initialization_method' = c('nndsvd', 'subSample'),
    'convergence_errors_print' = FALSE)
  
  eics_peaks.df <- as.data.frame(eics_peaks.mat)
  
  # test the main function
  total_res <- dia.nmf.f( rawData.onDiskMSnExp = rawData.onDiskMSnExp, ms1_peaks.df = eics_peaks.df,
             ppm = 7, rt_index = TRUE, 
             nmf_parameters.l = nmf_parameters.l, rank.method, initialization.method = "nndsvd",
             rt_tol = 0.1)
  
  
})