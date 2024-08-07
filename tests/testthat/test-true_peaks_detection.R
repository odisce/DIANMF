test_that("Extract the MS1 peaks", {
  if( !requireNamespace("xcms", quietly = TRUE) ){
    stop(
      "This package: xcms need to be installed to run the test."
    )
  }
  
  load("~/DIA_NMF_R_package/dianmf/data/data_example.rda")
  
  eics_rawData <- detect_EICs.f(
    rawData.onDiskMSnExp = data_example,
    d.out = NULL,
    ppm = 6,
    peakwidth = c(3,60), 
    snthresh = 0,
    prefilter = c(5,4000), 
    mzCenterFun = "wMeanApex3",
    integrate = 2,
    mzdiff = -0.001, 
    noise = 0, 
    firstBaselineCheck = FALSE
  )
  eics_peaks.mat <- xcms::chromPeaks(eics_rawData)
  
  expect_true("XCMSnExp" %in% class(eics_rawData))
  expect_true("matrix" %in% class(eics_peaks.mat))
  expect_true( nrow(eics_peaks.mat) >= 1 )
  
})
