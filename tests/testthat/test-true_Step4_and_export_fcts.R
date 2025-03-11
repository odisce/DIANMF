test_that("step 4 and export functions test", {
  
  features.l <- DIANMF.f(msexp = data_example,
                        d.out = FALSE,
                        sample_idx = 1, 
                        MS2_ISOEACHL = T, MS1MS2_L = F,
                        rank = 10,
                        maximumIteration = 200,
                        maxFBIteration = 100,
                        toleranceFB = 1e-05,
                        initialization_method = "nndsvd",
                        errors_print = FALSE,
                        method = "svsd",
                        scan_rt_ext = 10, min_distance = 5  )
  
  
  msSpectra <- exportMSSpectra(features.l[[1]])
    
  expect_true( all(sapply(msSpectra, function(x) inherits(x, "Spectra"))) )
  expect_true( !is.null(features.l) )
  
})