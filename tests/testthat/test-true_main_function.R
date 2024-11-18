test_that("Extract the MS1 peaks", {
  
  file <- system.file("extdata", "test_data.mzml", package = "DIANMF")
  ms1_peaks.df <- readRDS("~/DIA_NMF_R_package/dianmf/tests/testthat/testdata/ms1_peaks.rds")
  features.l_ms2 <- dia_nmf.f( mzML_path = file,
                           ms_level = "MS2",
                           peaks_by_xcms = FALSE, ms1_peaks = ms1_peaks.df,
                           ppm.n = 7,
                           maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                           MS1_init_method = 'nndsvd', MS2_init_method = 'subSample',
                           rt_tol = 2 )
  
  
  features.l_ms1 <- dia_nmf.f( mzML_path = file,
                           ms_level = "MS2",
                           peaks_by_xcms = TRUE, ms1_peaks = NULL,
                           ppm.n = 7,
                           maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                           MS1_init_method = 'nndsvd', MS2_init_method = 'subSample',
                           rt_tol = 2,
                           ppm = 6, peakwidth = c(3,60), snthresh = 0,
                           prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                           integrate = 2, mzdiff = -0.001, noise = 0,
                           firstBaselineCheck = FALSE)
  
  expect_equal(class(features.l_ms2), "list")
  expect_equal(length(features.l_ms2), 2)
  expect_true( !is.null(features.l_ms2[[1]]) ) 

  expect_equal(class(features.l_ms1), "list")
  expect_true( !is.null(features.l_ms1[[1]]) ) 

  
  p <- plot_ms2spectra_detailes(features.l_ms1[[1]]$MS2_pure_spectrum, features.l_ms1[[1]]$MS2_pure_spectrum_specific )
  
  # # I need to plot 3 spectra of Dextromorphan: NMF MS2 pure spectrum, dia spectrum and dda spectrum
  # 
  # nmf_spectrum_specific <- features.l_ms2[[1]]$MS2_pure_spectrum_specific
  # nmf_spectrum <- features.l_ms2[[1]]$MS2_pure_spectrum
  # clean_dia_database_5ppm <- readRDS("~/1workflow/DIA_NMF_workflow/Results/original_vs_clean_dia_databse/ppm_5/clean_dia_database_5ppm.rds")
  # dia_spectrum <- clean_dia_database_5ppm[[18]]
  # dda_spectrum <- sub_database$spectra[[1]]@spectrum
  # 
  # plot_spectra_vs(measured_spectrum = nmf_spectrum_specific, library_spectrum = dia_spectrum)
  # plot_spectra_vs(measured_spectrum = nmf_spectrum, library_spectrum = dia_spectrum)
  # plot_spectra_vs(measured_spectrum = nmf_spectrum, library_spectrum = dda_spectrum)
  
})