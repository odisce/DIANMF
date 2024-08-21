test_that("Extract the MS1 peaks", {
  if( !requireNamespace("xcms", quietly = TRUE) ){
    stop(
      "This package: xcms need to be installed to run the test."
    )
  }
  if( !requireNamespace("MSnbase", quietly = TRUE) ){
    stop(
      "This package: xcms need to be installed to run the test."
    )
  }
  
  load("~/DIA_NMF_R_package/dianmf/data/data_example.rda")
  ms1_peaks.mat <- detect_peaks_by_xcms(rawData.onDiskMSnExp = data_example,
                                         ppm = 6, peakwidth = c(3,60), snthresh = 0,
                                         prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                                         integrate = 2, mzdiff = -0.001, noise = 0,
                                         firstBaselineCheck = FALSE )
  
  ms1_peaks.df <- prepare_ms1_peaks(ms1_peaks = ms1_peaks.mat)
  
  expect_true("OnDiskMSnExp" %in% class(data_example))
  expect_true("matrix" %in% class(ms1_peaks.mat))
  expect_true("data.frame" %in% class(ms1_peaks.df))
  expect_true( nrow(ms1_peaks.df) >= 1 )
  
})
