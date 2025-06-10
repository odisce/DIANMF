test_that("Step 1 tests and some test for Step 2 and utils_function.R", {
 
  # test step 1-----------------------------------------------------------------
  input_dir <- system.file("extdata", package = "DIANMF")
  mzml_dt <- prepare_mzMLfiles(input_dir)
  expect_named( mzml_dt, c("mzml_path", "file_name", "creation_date", "class", "InjectionOrder") )

  mzml_seq <- create_seq(mzml_dt)
  expect_equal( nrow(mzml_dt), nrow(mzml_seq) )
  expect_named( mzml_seq, c("mzml_path", "class", "InjectionOrder") )

  params_ls <- list(
    "CentWaveParam" = xcms::CentWaveParam(
      ppm = 5,
      peakwidth = c(3, 10),
      snthresh = 2,
      prefilter = c(5, 1000),
      mzCenterFun = "wMeanApex3",
      integrate = 2,
      mzdiff = -0.005,
      noise = 100,
      firstBaselineCheck = FALSE
    ),
    "MergeNeighboringPeaksParam" = xcms::MergeNeighboringPeaksParam(
      expandRt = 2,
      expandMz = 0.001,
      ppm = 2,
      minProp = 0.75
    ),
    "ObiwarpParam" = xcms::ObiwarpParam(
      binSize = 0.05
    ),
    "PeakDensityParam" = xcms::PeakDensityParam(
      sampleGroups = NA,
      bw = 10,
      minFraction = 0.1,
      minSamples = 2,
      binSize = 0.01,
      ppm = 7,
      maxFeatures = 500
    ),
    "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
  )
  register(SnowParam(workers = 1))
  msexp <- detect_xcms_peaks(sequence_table = mzml_dt, params = params_ls)
  
  expect_true( 'XcmsExperiment' %in% class(msexp) )
  
  ms1_peaks <- extract_xcms_peaks(msexp = msexp)
  ms1_features_1 <- extract_xcms_features(msexp = msexp, orderL_sample = "test_data_20170803_FS-DIA-E2-10ng-rep3_pos_51.mzml", quantifyL = FALSE)
  ms1_features_2 <- extract_xcms_features(msexp = msexp, orderL_sample = NULL, quantifyL = TRUE)

  expect_true( "data.table" %in% class(ms1_peaks) )
  expect_true( "data.table" %in% class(ms1_features_1) & "data.table" %in% class(ms1_features_2) )
  expect_true( nrow(ms1_peaks) > 0 & nrow(ms1_features_1) > 0 & nrow(ms1_features_2) > 0 & 
                 nrow(ms1_features_1) == nrow(ms1_features_2) ) 


  # some utils_function.R tests ------------------------------------------------
  nev <- get_ms1_rtdiff(msexp) * 1.5
  expect_true( is.numeric(nev) )
  expect_true( nev > 0 )
  expect_true( all( c(1,2) %in% get_mslevels(msexp)) )
  
  raw_dt_i <- get_spectra_values(msexp)
  expect_true( "data.table" %in% class(raw_dt_i) )
  expect_true( nrow(raw_dt_i) > 0 )
  
  spect_info <- get_spectra_index(msexp)
  expect_true( nrow(spect_info) > 0 )
  
  # some tests from step 2-----------------------------------------------------------------
  
  ref <- 1
  ppm.n <- 10
  ref_range <- PpmRange(ref, ppm.n)
  expect_true(ref %between% ref_range )
  
  expect_true(has_n_consecutive_non_zero(c(0, 1, 1, 1, 1, 0), nscans = 4))
  expect_false(has_n_consecutive_non_zero(c(0, 1, 1, 0, 1, 0), nscans = 3))

  # ms1_features_peaks <- ms1_features_2[, .(peakindex = unlist(peakidx), mzmed, rtmed), by = featureid]
  # ms1_features_peaks <- merge(ms1_peaks, ms1_features_peaks, by = "peakindex")
  # 
  # expect_true(all(c("peakindex", "peakid", "sample", "is_filled", "featureid") %in% colnames(ms1_features_peaks)))
  # expect_true(nrow(ms1_features_peaks) > 0)
  # 
  # s_idx = 1
  # ms1_features <- copy(ms1_features_peaks[sample == s_idx, ])
  # ms1_features[, iteration := as.character(NA) ]
  # 
  # expect_true( unique(ms1_features$sample) == s_idx )
  # expect_true( all(is.na(ms1_features$iteration)) )
  # 
  # msexp_idx <- xcms::filterFile(msexp, s_idx)
  # file_info <- MsExperiment::sampleData(msexp_idx) %>% as.data.table()
  # spect_info <- get_spectra_index(msexp_idx)
  # 
  # expect_true( file_info$InjectionOrder == s_idx )
  # expect_true( class(msexp_idx) == "XcmsExperiment" )
  # expect_true( nrow(spect_info) > 0 )
  # expect_true( all( c(1,2) %in% get_mslevels(msexp_idx)) )
  # 
  # nev <- get_ms1_rtdiff(msexp_idx) * 1.5
  # 
  # expect_true( is.numeric(nev) )
  # expect_true( nev > 0 )
  # 
  # min_rt <- ms1_features[, min(rtmin)]
  # max_rt <- ms1_features[, max(rtmax)]
  # feature_idx <- 1
  # k <- 0
  # ## Sort and Extract features
  # ms1_features <- ms1_features[order(-into), ]
  # last_feature <- ms1_features[, last(featureid)]
  # 
  # feature_idx <- ms1_features[is_filled == 0 & is.na(iteration), ][1, featureid]
  # # feature_idx <- "FT04898"
  # peak_i <- ms1_features[featureid == feature_idx & is_filled == 0, ][which.max(into), ]
  # 
  # scan_rt_ext <- 10
  # # rt_method == "peak"
  # rt_range <- peak_i[, range(c(rtmin, rtmax)) + c(-scan_rt_ext - nev, +scan_rt_ext + nev)]
  # 
  # expect_true( is.numeric(rt_range) )
  # 
  # # rt_method == "constant"
  # rt_range <- peak_i[, rt + c(-scan_rt_ext - nev, +scan_rt_ext + nev)]
  # rt_range[1] <- max(min_rt, rt_range[1])  # to avoid rt < 0
  # rt_range[2] <- min(max_rt, rt_range[2])  # to avoid rt out of range
  # 
  # expect_true( is.numeric(rt_range) )
  # 
  # ## Subset msexp object
  # msexp_idx_rt <- xcms::filterRt(msexp_idx, rt_range) %>% suppressMessages()
  # 
  # rt_values <- range(rtime(msexp_idx_rt))
  # is_within <- rt_values[1] >= rt_range[1] && rt_values[2] <= rt_range[2]
  # expect_equal( is_within, TRUE )
  # 
  # ms1_peaks_i <- ms1_features[rtmin <= rt_range[2] & rtmax >= rt_range[1], -c("iteration")]
  # # flag peaks, and exclude feature not in rt range:
  # min_distance <- 3
  # border_lim <- min(c(min_distance, floor(diff(rt_range) / 3)))
  # rt_limits <- rt_range + c(+border_lim, -border_lim)
  # ms1_peaks_i[, peakfull := ifelse((rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
  #                                  ifelse(!rt %between% rt_limits, "apex_border",
  #                                         ifelse(rt %between% rt_range, "apex", "partial"
  #                                         )))]
  # ms1_peaks_i[, msLevel := 1]
  # ms1_peaks_i[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
  # ms1_peaks_i[, intensity := into]
  # features_iter <- ms1_peaks_i[peakfull %in% c("full", "apex"), unique(featureid)]
  # ## Exclude features when rtmed not in range
  # feat_notinR <- ms1_peaks_i[!rtmed %between% rt_limits, unique(featureid)]
  # features_iter <- features_iter %>% {.[!. %in% feat_notinR]} %>% c(., feature_idx) %>% unique()
  # ## Set all those features to '0' if not already extracted
  # ms1_features[featureid %in% features_iter & is.na(iteration), iteration := "0"]
  # 
  # # align scans and get time dictionary
  # time_dic <- align_scans(
  #   msexp = msexp_idx_rt,
  #   rt_range = rt_range,
  #   sample_idx = 1
  # )
  # 
  # expect_equal( colnames(time_dic), c("scan_norm", "isolationWindowTargetMz", "rtime", "msLevel") )
  # expect_true( 1 %in% unique(time_dic$msLevel) & 2 %in% unique(time_dic$msLevel) )
  # 
  # # normalized time dictionary
  # time_dic_i <- time_dic[, .(rtime = mean(rtime)), by = .(scan_norm)]
  # 
  # expect_equal( colnames(time_dic_i), c("scan_norm", "rtime") )
  # 
  # # Generating MS2 peaks list
  # combineSpectra_arg <- list(peaks = "intersect", ppm = 5, tolerance = 0.005)
  # ms2_peaks_i <- generate_peaklist(
  #   msexp = msexp_idx_rt,
  #   mslevel = 2,
  #   ms2isowinL = TRUE,
  #   combineSpectra_arg = combineSpectra_arg,
  #   ppm = 6
  # )
  # 
  # expect_equal( colnames(ms2_peaks_i), c("IsoWin", "msLevel", "mz", "rtime", "intensity",
  #                                        "isolationWindowTargetMz", "xic_label", "mzmin", "mzmax") )                  
  # expect_equal( unique(ms2_peaks_i$msLevel), 2 )
  # 
  # #  extract ions signals
  # raw_dt_i <- get_spectra_values(msexp_idx_rt)
  # 
  # expect_true( "data.table" %in% class(raw_dt_i) )
  # expect_true( nrow(raw_dt_i) > 0 )
  # 
  # # combine MS1 MS2 peaklists
  # MS1MS2_L <- TRUE
  # peaks_ls <- list(
  #   "MS1MS2" = rbind(
  #     ms1_peaks_i,
  #     ms2_peaks_i,
  #     fill = TRUE
  #   )
  # )
  # 
  # expect_equal( length(peaks_ls), 1 )
  # 
  # # Generating XICs
  # i <- length(peaks_ls)
  # ms_xics_i <- build_XICs(
  #   peaks_dt = peaks_ls[[i]],
  #   rawdt = raw_dt_i
  # )
  # ms_xics_i <- merge(ms_xics_i, time_dic[, .(rtime, scan_norm)], by = "rtime")
  # 
  # expect_true( nrow(ms_xics_i) > 0 & "data.table" %in% class(ms_xics_i) )
  # expect_equal( c("rtime", "msLevel", "isolationWindowTargetMz", "collisionEnergy", "xic_label",
  #                 "mz", "intensity", "scan_norm"  ), colnames(ms_xics_i) )
  # 
  # # Build mixed matrix
  # nscans <- 3
  # ms_mixed_i <- build_mixed_matrix(ms_xics_i, nscans = nscans)
  # 
  # expect_true( "matrix" %in% class(ms_mixed_i[[1]]) & nrow(ms_mixed_i[[1]]) > 0  )

})
