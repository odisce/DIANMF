test_that("nGMCAs algorithm test", {
  
  m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3)
  options.l <-  list(
    'maximumIteration' = 10,
    'maxFBIteration' = 5,
    'toleranceFB' = 1e-5,
    'useTranspose' = TRUE)
  ngmcas_res <- nGMCAs(originalMatrix = m, rank = 2, options = options.l, errors_print = TRUE, initialization_method = 'nndsvd')
  
  W <- ngmcas_res$S
  H <- ngmcas_res$A
  
  expect_true(all(W >= 0))
  expect_true(all(H >= 0))
  expect_equal(dim(W %*% H), dim(m))
  
})