test_that("Factorize mixed MS1 and MS2 data for a specific peak to get the pure spectra, 2 amino peak", {
  
  load("~/DIA_NMF_R_package/dianmf/data/data_example.rda")
  path_test <- testthat::test_path("testdata")
  features.l <- readRDS(paste0(path_test, "/features.rds"))
  
  ms1_mat <- features.l[[2]]$MS1_mixed_mat
  ms2_mat <- do.call(rbind, features.l[[2]]$MS2_mixed_mat)
  
  # MS1 data
  ngmcas_res <- nGMCAs(X.m = ms1_mat, rank = 3,
                       maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                       initialization_method = "nndsvd",
                       errors_print = FALSE)
  W_ms1 <- ngmcas_res$S
  H_ms1 <- ngmcas_res$A
  W_ms1 <- apply(W_ms1, 2, function(x) x / max(x))
  
  mz_prec <- features.l[[2]]$peak$mz
  rt_prec <- features.l[[2]]$peak$rt
  
  ms1_pure_data <- extract_ms1_pure_spectrum(W_ms1 = W_ms1, mz_prec = mz_prec)
  ms1_pure_spectrum <- ms1_pure_data$ms1_pure_spectrum
  comp_ms1 <- ms1_pure_data$comp_ms1

  # plot
  p_eics <- plot_MS_eics(ms_mixed = ms1_mat, ms_pure_H = H_ms1, ms_level = "MS1", rt_prec = rt_prec, choosen_comp = comp_ms1)
  p_spectra <- plot_MS_spectra(ms_mixed = ms1_mat, ms_pure_W = W_ms1, ms_level = "MS1", mz_prec, choosen_comp = comp_ms1)
  
  # MS2 data
  info.swath <- isolationWindows.range(data_example)
  
  rownames(H_ms1) <- NULL
  ngmcas_res_all <- nGMCAs(X.m = ms2_mat, rank = 3,
                           maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                           initialization_method = "subSample", H_sub = H_ms1,
                           errors_print = TRUE)
  W_ms2 <- ngmcas_res_all$S
  H_ms2 <- ngmcas_res_all$A
  W_ms2 <- apply(W_ms2, 2, function(x) x / max(x))
  
  comp_ms2 <- elutions_corelation(chromo_main = H_ms1[comp_ms1, ], chromos = H_ms2)
  ms2_pure_spectrum <- choose_ms2_pure_spectrum(W_ms2 = W_ms2, choosen_comp = comp_ms2)
  
  ms2_pure_spectrum_specific <- filter_ms2_spectrum(ms2_pure_spectrum = ms2_pure_spectrum,
                                                    ms2_matrices = features.l[[2]]$MS2_mixed_mat,
                                                    mz_prec = mz_prec, info.swath = info.swath )
  if( !is.null(ms2_pure_spectrum_specific) ){
    ms2_pure_spectrum_specific <- ms2_pure_spectrum_specific[ms2_pure_spectrum_specific['intensity'] != 0, ]
  }
  
  ms2_pure_spectrum <- ms2_pure_spectrum[ms2_pure_spectrum['intensity'] != 0, ]
  
  # plot
  p_eics <- plot_MS_eics(ms_mixed = ms2_mat, ms_pure_H = H_ms2, ms_level = "MS2", rt_prec = rt_prec, choosen_comp = comp_ms2)
  p_spectra <- plot_MS_spectra(ms_mixed = ms2_mat, ms_pure_W = W_ms2, ms_level = "MS2", mz_prec, choosen_comp = comp_ms2)
  
  
  # tests
  expect_true( !is.null(ms1_pure_spectrum) )
  expect_true( !is.null(ms2_pure_spectrum) )
  expect_true( !is.null(ms2_pure_spectrum_specific) )
  expect_equal(nrow(W_ms1), nrow(features.l[[2]]$W_ms1))  # tricky tests
  expect_equal(nrow(W_ms2), nrow(features.l[[2]]$W_ms2))
  expect_equal(ncol(H_ms1), ncol(features.l[[2]]$H_ms1))
  expect_equal(ncol(H_ms2), ncol(features.l[[2]]$H_ms2))
  
})