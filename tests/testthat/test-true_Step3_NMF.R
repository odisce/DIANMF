test_that("initializing methods test", {
  
  m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3)
  r <- 2
  
  # test random initialization
  random_res <- random_init(X = m, rank = r)
  W_random <- random_res$W
  H_random <- random_res$H
  
  expect_equal(dim(W_random %*% H_random), dim(m))
  expect_true(all(W_random >= 0))
  expect_true(all(H_random >= 0))
  expect_equal( ncol(W_random), r)
  expect_equal( nrow(H_random), r)
  
  
  # test nndsvd initialization
  nndsvd_res <- nndsvd_init(X = m, rank = r)
  W_nndsvd <- nndsvd_res$W
  H_nndsvd <- nndsvd_res$H
  
  expect_equal(dim(W_nndsvd %*% H_nndsvd), dim(m))
  expect_true(all(W_nndsvd >= 0))
  expect_true(all(H_nndsvd >= 0))
  expect_equal( ncol(W_nndsvd), r)
  expect_equal( nrow(H_nndsvd), r)
  
  
  # test sub-sampling initialization
  H_sub <- matrix(c(0.2,0,0.8,1,0.06,0.1), nrow = r)
  subSample_res <- subsample_init(Y = t(m), rank = r, H_sub = H_sub)
  A_subSample <- subSample_res$A
  S_subSample <- subSample_res$S
  
  expect_equal(dim( A_subSample %*% S_subSample), dim(t(m)))
  expect_true(all(S_subSample >= 0))
  expect_true(all(A_subSample >= 0))
  expect_equal( nrow(S_subSample), r)
  expect_equal( ncol(A_subSample), r)
  
  # test nGMCAs
  
  l1 <- get_svd_first(x = m, method = "base") 
  l2 <- get_svd_first(x = m, method = "fast") 
  l3 <- get_svd_first(x = m, method = "svds")
  expect_true( round(l1, 3) == round(l2, 3), round(l2, 3) == round(l3, 3))
  
  ngmcas_res <- nGMCAs(X.m = m, rank = 2,
                       maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                       initialization_method = 'nndsvd',
                       errors_print = FALSE, method = "fast")

  W <- ngmcas_res$S
  H <- ngmcas_res$A
  error <- error_function(Y = m, A = t(H), S = t(W))
  
  expect_true(all(W >= 0))
  expect_true(all(H >= 0))
  expect_equal(dim(W %*% H), dim(m))
  expect_true( error > 0 )
  
})
