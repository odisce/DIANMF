skip("development script")

test_that("test on barbier data 10ng replicate 3", {
  require(xcms)
  require(MSnbase)
  
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
  eics_peaks.mat <- xcms::chromPeaks(eics_peaks)
  
  expect_true("XCMSnExp" %in% class(eics_peaks))
  expect_true("matrix" %in% class(eics_peaks.mat))
  expect_true( nrow(eics_peaks.mat) >= 1 )
  
  eics_peaks.mat <- as.matrix(eics_peaks.mat[order(-eics_peaks.mat[, 'into']), ])
  idx.pg <-  1

  ms1_mat <- extract_ms1_matrix.f(idx.pg = idx.pg, eics_peaks.mat = eics_peaks.mat, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                  ppm = 7, rt_index = TRUE, mz_range = TRUE)
  
  expect_true("matrix" %in% class(ms1_mat))
  expect_true(dim(ms1_mat)[1] & dim(ms1_mat)[2] > 0)

  ms2_mat <- extract_ms2_matrix.f(idx.pg = idx.pg, eics_peaks.mat = eics_peaks.mat, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                  ppm = 7, rt_index = TRUE, mz_range = TRUE)

  expect_true("matrix" %in% class(ms2_mat))
  expect_true(dim(ms2_mat)[1] & dim(ms2_mat)[2] > 0)
  
  expect_equal(ncol(ms1_mat), ncol(ms2_mat))
  
  # test the nGMCAs algorithm
  
  
  
  
  
})
