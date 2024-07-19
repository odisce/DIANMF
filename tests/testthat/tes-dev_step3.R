skip("development script")

test_that("test on barbier data 10ng replicate 3; step3; extract mixed matrix for a specific peak", {
  devtools::load_all()
  
  file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
  rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");
  
  eics_peaks <- detect_EICs.f(
    rawData.onDiskMSnExp = rawData.onDiskMSnExp,
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
  expect_true("XCMSnExp" %in% class(eics_peaks))
  
  eics_peaks.mat <- xcms::chromPeaks(eics_peaks)
  eics_peaks.mat <- as.matrix(eics_peaks.mat[order(-eics_peaks.mat[, 'into']), ])
  idx.pg <-  1
  ms1_mat <- extract_ms1_matrix.f(idx.pg, eics_peaks.mat, rawData.onDiskMSnExp, ppm =15, rt_index = 1, mz_range = 1)
  
  expect_true("matrix" %in% class(ms1_mat))
  expect_true( nrow(ms1_mat) >= 1 )
})


# where I can add the R pcakage that are needed to run this script