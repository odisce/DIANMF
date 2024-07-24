test_that("initializing methods test", {
  
  m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3)
  r <- 2
  
  # test random initialization
  random_res <- random_init(A = m, k = r)
  W_random <- random_res$W
  H_random <- random_res$H
  
  expect_equal(dim(W_random %*% H_random), dim(m))
  expect_true(all(W_random >= 0))
  expect_true(all(H_random >= 0))
  expect_equal( ncol(W_random), r)
  expect_equal( nrow(H_random), r)
  
  
  # test nndsvd initialization
  nndsvd_res <- nndsvd_init(A = m, k = r)
  W_nndsvd <- nndsvd_res$W
  H_nndsvd <- nndsvd_res$H
  
  expect_equal(dim(W_nndsvd %*% H_nndsvd), dim(m))
  expect_true(all(W_nndsvd >= 0))
  expect_true(all(H_nndsvd >= 0))
  expect_equal( ncol(W_nndsvd), r)
  expect_equal( nrow(H_nndsvd), r)
  
  
  # test sub-sampling initialization
  H_sub <- matrix(c(0.2,0,0.8,1,0.06,0.1), nrow = r)         # I want to change the parameters of this function
  subSample_res <- subsample_init(mat = t(m), k = r, H_sub = t(H_sub))
  W_subSample <- subSample_res$W
  H_subSample <- subSample_res$H
  
  expect_equal(dim( H_subSample %*% W_subSample), dim(t(m)))
  expect_true(all(W_subSample >= 0))
  expect_true(all(H_subSample >= 0))
  expect_equal( nrow(W_subSample), r)
  expect_equal( ncol(H_subSample), r)
})
