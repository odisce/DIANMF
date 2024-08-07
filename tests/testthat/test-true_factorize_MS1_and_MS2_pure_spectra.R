test_that("Factorize mixed MS1 and MS2 data for a specific peak to get the pure spectra", {
  
  path_test <- testthat::test_path("testdata")
  features.l <- readRDS(paste0(path_test, "/features.rds"))
  
  MS1_mat <- features.l[[2]]$MS1_mixed_mat
  MS2_mat <- do.call(rbind, features.l[[2]]$MS2_mixed_mat)
  
  # Factorize MS1 data
  ngmcas_res <- nGMCAs(originalMatrix = MS1_mat, rank = 3,
                       maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                       initialization_method = "nndsvd",
                       errors_print = FALSE)
  W_ms1 <- ngmcas_res$S
  H_ms1 <- ngmcas_res$A
  
  # extract pure MS1 spectra and choose the peak corresponding spectrum
  mz_prec <- features.l[[2]]$peak$mz
  rt_prec <- features.l[[2]]$peak$rt
  
  mz_ms1_ions <- round(as.numeric(rownames(W_ms1)), 4);
  ms1_pure_spectra <- prepare_pure_spectra(W_ms1)
  
  closest_row <- which.min(abs(mz_ms1_ions - mz_prec));
  choosen_comp <- which.max(W_ms1[closest_row, ]);
  ms1_pure_spectrum <- ms1_pure_spectra[ms1_pure_spectra$comp_nb == paste0('comp',choosen_comp), ]
  ms1_pure_spectrum <- ms1_pure_spectrum[, c('mz_value', 'intensity')]
  ms1_pure_spectrum$intensity <- ms1_pure_spectrum$intensity / max(ms1_pure_spectrum$intensity)
  ms1_pure_spectrum <- ms1_pure_spectrum[ms1_pure_spectrum['intensity'] != 0, ]

  # factorize MS2 data
  rownames(H_ms1) <- NULL
  ngmcas_res_all <- nGMCAs(originalMatrix = MS2_mat, rank = 3,
                           maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                           initialization_method = "subSample", H_sub = H_ms1,
                           errors_print = FALSE)
  W_ms2 <- ngmcas_res_all$S
  H_ms2 <- ngmcas_res_all$A
  
  ms2_pure_spectra <- prepare_pure_spectra(W_ms2)
  ms2_pure_spectra <-  ms2_pure_spectra %>%
    group_by(comp_nb) %>%
    mutate(intensity = intensity / max(intensity)) 
  
  ms2_pure_spectrum <- subset(ms2_pure_spectra, comp_nb == paste0('comp',choosen_comp));
  ms2_pure_spectrum <- ms2_pure_spectrum[, c('mz_value','intensity')]
  ms2_pure_spectrum <- ms2_pure_spectrum[ms2_pure_spectrum['intensity'] != 0, ]

  expect_true( !is.null(ms1_pure_spectrum) )
  expect_true( !is.null(ms2_pure_spectrum) )
  
})