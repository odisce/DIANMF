test_that("step 4 test", {
  
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
  
  expect_true( !is.null(features.l) )
  
})