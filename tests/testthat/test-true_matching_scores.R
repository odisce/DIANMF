test_that("scores matching test", {
  
  spectrum1 <- data.frame(
    "mz_value" = seq(1:8),
    'intensity' = c(0,400,1000,200,0,700,300,300)
  )
  spectrum2 <- data.frame(
    "mz_value" = seq(1:8),
    'intensity' = c(100,0,1000,0,100,0,400,300)
  ) # refernce spectra
 
  dot_product <- GetSimpleDotProductSimilarity(measuredSpectra = spectrum1, librarySpectra = spectrum2, bin = 0.001)
  reverse_dot_product <- getReverseSearchingSimilarity(measuredSpectra = spectrum1, librarySpectra = spectrum2, bin = 0.001)
  presence_ratio <- GetPresenceSimilarity(measuredSpectra = spectrum1, librarySpectra = spectrum2, bin = 0.001)
  
  expect_true(is.numeric(dot_product))
  expect_equal( round(dot_product,3), 0.528 )
  
  expect_true(is.numeric(reverse_dot_product))
  expect_equal( round(reverse_dot_product,3), 0.872 )
  
  expect_true(is.numeric(presence_ratio))
  expect_equal( round(presence_ratio,3), 0.6 )
  
})
# spectrum1 <- data.frame(
#   "mz_value" = seq(1:5),
#   'intensity' = c(0,1000,0,300,300)
# )
# spectrum2 <- data.frame(
#   "mz_value" = seq(1:5),
#   'intensity' = c(100,1000,100,400,300)
# ) # refernce spectra
