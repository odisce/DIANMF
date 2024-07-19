test_that("load data test", {
  if( !requireNamespace("msdata", quietly = TRUE) ){
    stop(
      "This package: msdata need to be installed to run the test."
    )
  }
  if( !requireNamespace("xcms", quietly = TRUE) ){
    stop(
      "This package: xcms need to be installed to run the test."
    )
  }
  if( !requireNamespace("MSnbase", quietly = TRUE) ){
    stop(
      "This package: MSnbase need to be installed to run the test."
    )
  }
  
  mzml_path <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML", package = "msdata")
  rawData.onDiskMSnExp <- MSnbase::readMSData(mzml_path, mode = "onDisk");
  
  eics_rawData <- detect_EICs.f(
    rawData.onDiskMSnExp = rawData.onDiskMSnExp,
    d.out = NULL,
    ppm = 6,
    peakwidth = c(2,60), 
    snthresh = 0,
    prefilter = c(2,10), 
    mzCenterFun = "wMeanApex3",
    integrate = 2,
    mzdiff = -0.001, 
    noise = 0, 
    firstBaselineCheck = FALSE
  )
  eics_peaks <- xcms::chromPeaks(eics_rawData)
  
  expect_true("XCMSnExp" %in% class(eics_rawData))
  expect_true("matrix" %in% class(eics_peaks))
  expect_true( nrow(eics_peaks) >= 1 )
  
})
