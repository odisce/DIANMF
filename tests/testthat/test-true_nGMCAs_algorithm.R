test_that("nGMCAs algorithm test", {
  
  m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3)
  ngmcas_res <- nGMCAs(originalMatrix = m, rank = 2,
                       maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                       initialization_method = 'nndsvd',
                       errors_print = TRUE)
  
  W <- ngmcas_res$S
  H <- ngmcas_res$A
  
  expect_true(all(W >= 0))
  expect_true(all(H >= 0))
  expect_equal(dim(W %*% H), dim(m))
  
})