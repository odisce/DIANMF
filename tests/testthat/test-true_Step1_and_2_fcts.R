test_that("Step 1 and 2 tests", {
 
  # test step 1
  data_example <- get_test_peaks()
  expect_true( 'XcmsExperiment' %in% class(data_example) ) # test a data example function
  
  ms1_peaks <- extract_xcms_peaks(msexp = data_example)
  ms1_features <- extract_xcms_features(data_example)
  
  
  expect_true( isTRUE( "data.table" %in% class(ms1_peaks) ) )
  expect_true( isTRUE( "data.table" %in% class(ms1_features) )   )
  expect_true( nrow(ms1_features) > 0 & nrow(ms1_peaks) > 0 )
  
  # test step 2
  sample_idx <- 1
  rt_range <- c(537, 560)
  ms1_peaks <- extract_xcms_peaks(data_example)
  
  res_general <- get_rawD_ntime(msexp = data_example, rt_range, sample_idx)
  raw_dt <- res_general$raw_dt
  time_dic <- res_general$time_dic
  
  peaks_i <- ms1_peaks[sample == sample_idx & rtmin <= rt_range[2] & rtmax >= rt_range[1], ]
  peaks_i[, peakfull := ifelse(
    (rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
    ifelse(rt %between% rt_range, "apex", "partial")
  )]
  
  xic_dt_ms1 <- build_ms1XICS(peaks_i, raw_dt)
  xic_dt_ms2 <- build_ms2XICs(msexp = data_example, raw_dt, time_dic, rt_range, MS2_ISOEACHL = T)
  
  ### ms1
  ms1Data <- ms1Info(xic_dt_ms1)
  ms1_mixedmat <- ms1Data$ms1_mixedmat
  ms1_mixedmat_deleted <- ms1Data$ms1_mixedmat_deleted
  ms1_infos <- ms1Data$ms1_infos
  
  ### ms2
  ms2Data <- ms2Info(xic_dt_ms2)
  ms2_mixedmat <- ms2Data$ms2_mixedmat
  ms2_infos <- ms2Data$ms2_infos
  
  expect_true( "matrix" %in% class(ms2_mixedmat) &&  "matrix" %in% class(ms1_mixedmat) )
  expect_true( nrow(ms2_mixedmat) > 0 &&   nrow(ms1_mixedmat) > 0  )
  expect_equal( ncol(ms2_mixedmat), ncol(ms1_mixedmat) )
  
})
