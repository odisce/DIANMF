test_that("initializing methods test", {
  
  # test random initialization
  m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3)
  r <- 2
  random_res <- random_init(A = m, rank = r)
  W_random <- random_res$W
  H_random <- random_res$H
  
  expect_equal(dim(W_random %*% H_random), dim(m))
  expect_true(all(W_random >= 0))
  expect_true(all(H_random >= 0))
  expect_equal( ncol(W_random), r)
  expect_equal( nrow(H_random), r)
  
  
  # test nndsvd initialization
  nndsvd_res <- nndsvd(A = m, k = r)
  W_nndsvd <- nndsvd_res$W
  H_nndsvd <- nndsvd_res$H
  
  expect_equal(dim(W_nndsvd %*% H_nndsvd), dim(m))
  expect_true(all(W_nndsvd >= 0))
  expect_true(all(H_nndsvd >= 0))
  expect_equal( ncol(W_nndsvd), r)
  expect_equal( nrow(H_nndsvd), r)
  
  
  # test sub-sampling initialization
  # soon...
})
