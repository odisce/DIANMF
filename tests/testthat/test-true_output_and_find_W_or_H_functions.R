test_that("Given X and W, find H or the inverse", {
  
  path_test <- testthat::test_path("testdata")
  features.l <- readRDS(paste0(path_test, "/features.rds"))
  
  # test functions that extract information from the output object
  # extract pure MS1 spectra
  MS1_spect <- extract_pureSpect(features = features.l, spec_level = 'MS1')
  # extract pure MS2 spectra
  MS2_spect <- extract_pureSpect(features = features.l, spec_level = 'MS2')
  # extract pure MS1 spectra
  MS2_spectPrec2 <- extract_pureSpect(features = features.l, spec_level = 'MS2_specific')
  
  # ms1 raw matrix
  ms1_mat <- extract_mixedMat(features = features.l, ms_level = 'MS1')
  # ms2 raw matrix
  ms2_mat <- extract_mixedMat(features = features.l, ms_level = 'MS2')
  
  # factorized matrices
  H_ms1 <- extract_pureMat(features = features.l, ms_level = 'MS1', H = TRUE)
  W_ms1 <- extract_pureMat(features = features.l, ms_level = 'MS1', H = FALSE)
  H_ms2 <- extract_pureMat(features = features.l, ms_level = 'MS2', H = TRUE) 
  W_ms2 <- extract_pureMat(features = features.l, ms_level = 'MS2', H = FALSE)
  
  #-----------------------------------------------------------------------------
  
  # test functions to fund W or H
  H1 <- find_H(X = ms1_mat[[1]], W = W_ms1[[1]], maxFBIteration = 10, toleranceFB = 1e-5) # for ms1 data
  H2 <- find_H(X = ms2_mat[[1]], W = W_ms2[[1]], maxFBIteration = 10, toleranceFB = 1e-5) # for ms2 data
  
  W1 <- find_W(X = ms1_mat[[1]], H = H_ms1[[1]], maxFBIteration = 10, toleranceFB = 1e-5) # for ms1 data
  W2 <- find_W(X = ms2_mat[[1]], H = H_ms2[[1]], maxFBIteration = 10, toleranceFB = 1e-5) # for ms2 data
  
  expect_equal(dim( W1 %*% H1), dim(ms1_mat[[1]]))
  expect_equal(dim( W2 %*% H2), dim(ms2_mat[[1]]))
  
})