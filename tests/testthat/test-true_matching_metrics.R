test_that("matching metrics test", {
  
  spectrum1 <- data.frame(
    "mz_value" = seq(1:8),
    'intensity' = c(0,400,1000,200,0,700,300,300)
  ) # measured spectrum
  spectrum2 <- data.frame(
    "mz_value" = seq(1:8),
    'intensity' = c(100,0,1000,0,100,0,400,300)
  ) # reference spectrum
  
  tolerance <- 0.007
  
  dp_score <- GetSimpleDotProductSimilarity(measuredSpectra = spectrum1, librarySpectra = spectrum2, bin = tolerance)
  rdp_score <- getReverseSearchingSimilarity(measuredSpectra = spectrum1, librarySpectra = spectrum2, bin = tolerance)
  fpp_score <- GetPresenceSimilarity(measuredSpectra = spectrum1, librarySpectra = spectrum2, bin = tolerance)
  
  expect_true(is.numeric(dp_score))
  expect_equal( round(dp_score,3), 0.523 )
  
  expect_true(is.numeric(rdp_score))
  expect_equal( round(rdp_score,3), 0.561 )
  
  expect_true(is.numeric(fpp_score))
  expect_equal( round(fpp_score,3), 0.6 )
  
  
  # test the matching function with library database
  
  
  
})