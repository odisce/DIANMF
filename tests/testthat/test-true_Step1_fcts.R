test_that("Step 1 tests", {
 
  ms1_peaks <- extract_xcms_peaks(msexp = data_example)
  ms1_features <- extract_xcms_features(data_example)
  
  
  expect_true( isTRUE( "data.table" %in% class(ms1_peaks) ) )
  expect_true( isTRUE( "data.table" %in% class(ms1_features) )   )
  expect_true( nrow(ms1_features) > 0 & nrow(ms1_peaks) > 0 )
  
})