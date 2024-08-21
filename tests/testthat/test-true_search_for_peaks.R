test_that("scores matching test", {
  
  # I added the Barbier_supplementary_tableS1.tsv in inst/extdata 
  compounds_tb <- system.file("extdata", "Barbier_supplementary_tableS1.tsv", package = "DIANMF")
  ms1_peaks.df <- readRDS("~/DIA_NMF_R_package/dianmf/tests/testthat/testdata/ms1_peaks.rds")
  
  compounds.df <- load_compounds_function.tb(compounds_tb)
  
  peaks <- search_for_ROIs_function.l(compounds.df, ms1_peaks = ms1_peaks.df, mz_tol.n = 0.05, rt_tol.n = 5)

  expect_equal( nrow(compounds.df), 50 )
  expect_true(class(compounds.df) == "data.frame")
  
})
