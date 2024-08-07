skip("development script")

test_that("Create the unitary tests data", {
  
  load("~/DIA_NMF_R_package/dianmf/data/data_example.rda")
  eics_rawData <- detect_EICs.f(
    rawData.onDiskMSnExp = data_example,
    d.out = NULL,
    ppm = 6,
    peakwidth = c(3,60), 
    snthresh = 0,
    prefilter = c(5,4000), 
    mzCenterFun = "wMeanApex3",
    integrate = 2,
    mzdiff = -0.001, 
    noise = 0, 
    firstBaselineCheck = FALSE
  )
  eics_peaks.mat <- xcms::chromPeaks(eics_rawData)
  saveRDS(eics_peaks.mat, 'C:/Users/DK273056/Documents/DIA_NMF_R_package/dianmf/tests/testthat/testdata/eics_peaks_mat.rds')

  features.l <- dia.nmf.f( rawData.onDiskMSnExp = data_example, ms1_peaks.mat = eics_peaks.mat,
                           peaks_nb = 2,
                           ppm = 7, 
                           maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                           MS1_init_method = 'nndsvd', MS2_init_method = 'subSample', errors_print = FALSE,
                           rt_tol = 2)
  saveRDS(features.l, 'C:/Users/DK273056/Documents/DIA_NMF_R_package/dianmf/tests/testthat/testdata/features.rds')
  
})