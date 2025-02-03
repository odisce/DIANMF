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
  MSobject <- get_test_peaks()
  ms1_peaks.mat <- extract_xcms_peaks(
    MSobject,
    sample_nb = 1
  )
  
  ms1_peaks.df <- prepare_ms1_peaks(ms1_peaks = ms1_peaks.mat)
  
  expect_true("OnDiskMSnExp" %in% class(data_example))
  expect_true("matrix" %in% class(ms1_peaks.mat))
  expect_true("data.frame" %in% class(ms1_peaks.df))
  expect_true( nrow(ms1_peaks.df) >= 1)
})
