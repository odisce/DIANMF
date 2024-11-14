skip("development script")

test_that("test main on barbier data", {
  require(MSnbase)
  require(xcms)

  file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
  rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");
  fast_sample <- FALSE
  if (fast_sample) {
    sample_file <- MSnbase::filterRt(rawData.onDiskMSnExp, c(0,2*60))
  } else {
    sample_file <- rawData.onDiskMSnExp
  }
  eics_peaks.mat <- detect_peaks_by_xcms(rawData.onDiskMSnExp = sample_file,
                                     ppm = 7, peakwidth = c(6,60), snthresh = 1,
                                     prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                                     integrate = 2, mzdiff = -0.001, noise = 0,
                                     firstBaselineCheck = FALSE )
  
  # test the main function
  total_res <- dia_nmf.f( mzML_path = file,
                          ms_level = "MS1",
                          peaks_by_xcms = FALSE, ms1_peaks = eics_peaks.mat,
                          ppm.n = 7,
                          maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
                          MS1_init_method = 'nndsvd', MS2_init_method = 'subSample', errors_print = FALSE,
                          rt_tol = 2 )
  
  total_res2 <- dia_nmf.f(mzML_path = file,
                          ms_level = "MS2",
                          peaks_by_xcms = FALSE, ms1_peaks = eics_peaks.mat[1:35, ], d.out = NULL,
                          ppm.n = 7,
                          maximumIteration = 5, maxFBIteration = 10, toleranceFB = 1e-5,
                          MS1_init_method = 'nndsvd', MS2_init_method = 'subSample', errors_print = FALSE,
                          rt_tol = 2,
                          ppm = 7, peakwidth = c(6,60), snthresh = 1,
                          prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
                          integrate = 2, mzdiff = -0.001, noise = 0,
                          firstBaselineCheck = FALSE)
  
  # eics_peaks.mat[order(-eics_peaks.mat[, "into"]),][1:10, ]
  # plot_path = "C:/Users/DK273056/Documents/DIA_NMF_R_package_plots"
  # saveRDS(total_res2, 'C:/Users/DK273056/Documents/DIA_NMF_R_package_outputs/features.rds')
  
 # # to analyse the running time
 # running_time <- profvis::profvis({
 #   res_sub <- dia.nmf.f(  )
 #   },
 #   interval = 0.1
 #   )
 # 
 # # Save a profile to an HTML file
 # htmlwidgets::saveWidget(running_time, "profile.html")
 # # Can open in browser from R
 # browseURL("profile.html")
 # [order(-eics_peaks.mat[, "into"]),][1:10, ]
})