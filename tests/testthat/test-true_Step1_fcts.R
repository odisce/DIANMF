test_that("Step 1 tests", {
 
  seq <- create_seq("./inst/extdata/test_data1.mzml")
  ms1_peaks <- extract_xcms_peaks(msexp = data_example)
  ms1_features <- extract_xcms_features(data_example)
  
  
  expect_true( ncol(seq) == 3 )
  expect_equal( colnames(seq), c("mzml_path", "class" ,"InjectionOrder") )
  expect_true( isTRUE( "data.table" %in% class(ms1_peaks) ) )
  expect_true( isTRUE( "data.table" %in% class(ms1_features) )   )
  expect_true( nrow(ms1_features) > 0 & nrow(ms1_peaks) > 0 )
  
})
