skip("development script")

test_that("Create the unitary tests data", {
  
  load("~/DIA_NMF_R_package/dianmf/data/data_example.rda")
  ms1_peaks.mat <- detect_peaks_by_xcms(rawData.onDiskMSnExp = data_example,
                                        ppm = 6, peakwidth = c(3,60), snthresh = 0,
                                        prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                                        integrate = 2, mzdiff = -0.001, noise = 0,
                                        firstBaselineCheck = FALSE )
  
  ms1_peaks.df <- prepare_ms1_peaks(ms1_peaks = ms1_peaks.mat)
  
  saveRDS(ms1_peaks.df[1:2, ], 'C:/Users/DK273056/Documents/DIA_NMF_R_package/dianmf/tests/testthat/testdata/ms1_peaks.rds')
  
  file <- system.file("extdata", "test_data.mzml", package = "DIANMF")
  features.l <- dia_nmf.f( mzML_path = file,
                           ms_level = "MS2",
                           peaks_by_xcms = FALSE, ms1_peaks = ms1_peaks.df[1:2, ],
                           ppm.n = 7,
                           maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                           MS1_init_method = 'nndsvd', MS2_init_method = 'subSample',
                           rt_tol = 2 )
  saveRDS(features.l, 'C:/Users/DK273056/Documents/DIA_NMF_R_package/dianmf/tests/testthat/testdata/features.rds')
  
  # features.l2 <- dia_nmf.f( mzML_path = file,
  #                          ms_level = "MS1",
  #                          peaks_by_xcms = TRUE, ms1_peaks = NULL, d.out = NULL,
  #                          ppm = 6, peakwidth = c(3,60), snthresh = 0,
  #                          prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
  #                          integrate = 2, mzdiff = -0.001, noise = 0,
  #                          firstBaselineCheck = FALSE)
  
})