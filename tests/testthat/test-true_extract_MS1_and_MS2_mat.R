test_that("Extract mixed MS1 and MS2 data for a specific peak", {
  
  load("~/DIA_NMF_R_package/dianmf/data/data_example.rda")
  
  path_test <- testthat::test_path("testdata")
  eics_peaks.mat <- readRDS(paste0(path_test, "/eics_peaks_mat.rds"))
  
  eics_peaks.mat <- as.matrix(eics_peaks.mat[order(-eics_peaks.mat[, 'into']), ])
  ms1_peaks.df <- as.data.frame(eics_peaks.mat)
  
  # this peak is related to Aminophenol spiked precursor------------------------
  peak.idx <- 2 
  
  # extract the MS1 mixed data--------------------------------------------------
  ms1_mat <- extract_ms_matrix.f(idx.pg = peak.idx, eics_peaks.mat = ms1_peaks.df,
                                 rawData.onDiskMSnExp = data_example,
                                 ppm = 7, rt_index = TRUE, mz_range = NULL, iso_win_index = NULL)
  row_filter <- apply(ms1_mat, 1, has_four_consecutive_non_zero)
  ms1_mat <- ms1_mat[row_filter, , drop = FALSE]
  
  expect_true( nrow(ms1_mat) >= 1 )
  expect_true( ncol(ms1_mat) > 1 )

  # extract the MS2 data--------------------------------------------------------
  idx.swath <- seq(1,3)
  ms2.l <- lapply(idx.swath, function(i){
    ms2_mat <- extract_ms_matrix.f(idx.pg = peak.idx, eics_peaks.mat = ms1_peaks.df,
                                   rawData.onDiskMSnExp = data_example,
                                   ppm = 7, rt_index = TRUE, mz_range = NULL, iso_win_index = i)
    return(ms2_mat)
  })
  ms2.l <- map(ms2.l, filter_matrix_rows)
  ms2_matrices <- do.call(rbind, ms2.l)
  
  expect_true( nrow(ms2_matrices) >= 1 )
  expect_true( ncol(ms2_matrices) > 1 )
  
  expect_equal(ncol(ms1_mat), ncol(ms2_matrices))

})
