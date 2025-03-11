testthat::skip(" to be solved soon  ")
test_that("data example functions test", {
  library(stats)

  mzml_files_paths <- get_test_mzml()
  mzml_dt_data_test <- get_test_sequence()
  # xcms_obj <- get_test_peaks()

  expect_equal(length(mzml_files_paths), 2)
  expect_true( "data.frame" %in% class(mzml_dt_data_test) )
  expect_true( ncol(mzml_dt_data_test) == 3 )
  expect_equal( colnames(mzml_dt_data_test), c("mzml_path", "class" ,"InjectionOrder") )
  # expect_true( 'XcmsExperiment' %in% class(xcms_obj) )

})
