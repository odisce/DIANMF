test_that("data example functions test", {

  mzml_files_paths <- get_test_mzml()
  mzml_dt_data_test <- get_test_sequence()
  mzml_peaks <- get_test_peaks()
  
  expect_equal(length(mzml_files_paths), 2)
  expect_true( "data.frame" %in% class(mzml_dt_data_test) )
  expect_true( ncol(mzml_dt_data_test) == 3 )
  expect_named( mzml_dt_data_test, c("mzml_path", "class" ,"InjectionOrder") )
  expect_true( "XcmsExperiment" %in%  class(mzml_peaks) )

  file_info <- MsExperiment::sampleData(mzml_peaks) %>% as.data.table()
  expect_equal( nrow(file_info), 2 )
  expect_true( "A" %in% file_info$class )
  expect_true( all(file_info$InjectionOrder == c(1,2)) )

})
