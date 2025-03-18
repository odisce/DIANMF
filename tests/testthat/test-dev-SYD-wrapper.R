testthat::skip("development script; test the wrapper function")

{
  devtools::document()
  devtools::load_all()
  input_dir <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
  input_targ <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Compound_data_corrected_rt.csv"
  targ_dt <- fread(input_targ)
  mzml_dt <- prepare_mzMLfiles(input_dir = input_dir)
  mzml_dt <- mzml_dt[class == 10.00, ]
  params_ls <- list(
    "CentWaveParam" = xcms::CentWaveParam(
      ppm = 6,
      peakwidth = c(6, 30),
      snthresh = 2,
      prefilter = c(5, 4000),
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
      binSize = 0.008,
      ppm = 7,
      maxFeatures = 500
    ),
    "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
  )

  cache_path <- "/spi/scidospi/06_Data/deconv_nmf/temp.rds"
  profile_path <- "/spi/scidospi/06_Data/deconv_nmf"
  if (file.exists(cache_path)) {
    xcms_obj <- readRDS(cache_path)
  } else {
    xcms_obj <- DIANMF::detect_xcms_peaks(sequence_table = mzml_dt, params = params_ls)
    saveRDS(xcms_obj, cache_path)
  }
}

{
  require(profvis)
  devtools::load_all()
  profile_save_path <- file.path(profile_path, format(Sys.time(), "%y%m%d-%H%M-profvis.html"))
  proffres <- profvis::profvis(
    expr = {
      features <- DIANMF.f(
        msexp = xcms_obj, dir_out = FALSE,
        sample_idx = 1,
        MS2_ISOEACHL = TRUE,
        MS1MS2_L = FALSE,
        rank = 10,
        maximumIteration = 200,
        maxFBIteration = 100,
        toleranceFB = 1e-05,
        initialization_method = "nndsvd",
        errors_print = FALSE,
        method = "svds",
        scan_rt_ext = 10,
        min_distance = 5,
        featuresn = 2
      )
    }
  )
  htmlwidgets::saveWidget(proffres, profile_save_path)
}

devtools::load_all()
{
  cache_path <- "/spi/scidospi/06_Data/deconv_nmf/features.rds"
  if (file.exists(cache_path)) {
    features <- readRDS(cache_path)
  } else {
    data.table::setDTthreads(1)
    A <- Sys.time()
    rt_range <- targ_dt[grepl("Scopolamine", Compound), (rt_sec + c(-60, +60))]
    msexp <- xcms::filterRt(xcms_obj, rt_range)
    # msexp <- xcms_obj
    features <- DIANMF.f(
      msexp = msexp,
      dir_out = FALSE,
      sample_idx = NULL,
      MS2_ISOEACHL = TRUE,
      MS1MS2_L = FALSE,
      rank = 20,
      min_contrib = 0.6,
      maximumIteration = 200,
      maxFBIteration = 100,
      toleranceFB = 1e-05,
      initialization_method = "nndsvd",
      errors_print = FALSE,
      method = "svds",
      scan_rt_ext = 15,
      min_distance = 4,
      featuresn = NULL,
      nscans = 6,
      rt_method = "constant",
      combineSpectra_arg = list(peaks = "intersect", ppm = 4, tolerance = 0.005, minProp = 0.05),
      verbose = TRUE
    )
    B <- Sys.time()
    difftime(B, A)
    # 2.28 hours for 3 files
    saveRDS(features, cache_path)
  }
}

{
  devtools::load_all()
  ## Features summary
  ### Annotate features
  temp_ft <- get_feature_summary(features.l = features)
  ft_match <- search_features(
    feature_dt = temp_ft,
    dt = targ_dt[, .(Compound, mz = mz_pos, rt = rt_sec)],
    rttol = 20,
    ppm = 8
  )
  temp_ft <- temp_ft[order(-into), ]
  ## Plot one sample diagnostic plot for one iteration
  save_path <- "/spi/scidospi/06_Data/deconv_nmf/20250318_test"
  dir.create(save_path)
  # list.files(save_path, full.names = TRUE) %>% sapply(., unlink)
  for (feati in ft_match[, unique(featureid)]) {
    timeA <- Sys.time()
    message(sprintf("Extracting: %s", feati), appendLF = FALSE)
    plot_out <- lapply(temp_ft[, unique(sample)], function(spli) {
      plot_feature(
        features.l = features,
        summary_dt = temp_ft,
        feature_id = feati,
        sample_index = spli,
        log2L = FALSE,
        max_method = "contribution"
      ) %>%
        suppressWarnings()
    }) %>% {
      ggpubr::ggarrange(plotlist = ., ncol = length(.))
    }
    ggsave(
      sprintf(
        "%s/%s.png",
        save_path,
        feati
      ),
      plot = plot_out,
      bg = "white",
      w = 16,
      h = 9
    )
    timeB <- Sys.time()
    message(
      sprintf(
        ": Done in %0.2fs",
        difftime(timeB, timeA, units = "secs")
      ),
      appendLF = TRUE
    )
  }
  
  
  ## Plot pure MS1/MS2 for a feature in different samples (with MS1 & MS2 matching)

  ## Plot xcms peaks ranges
  require(ggplot2)
  ggplot(ms1_features[is_filled == 0,]) +
    geom_linerange(aes(y = mz, x = rt, xmin = rtmin, xmax = rtmax)) +
    geom_linerange(data = ms1_peaks_i[is_filled == 0,], aes(y = mz, x = rt, xmin = rtmin, xmax = rtmax, color = peakfull)) + 
    theme_bw()

  nmf_result_ls$pure_rt[, value_sc := scales::rescale(value, to = c(0,1)), by = .(MSid, rank)]
  ggplot(nmf_result_ls$pure_rt, aes(rtime, value_sc, color = rank, group = rank)) +
    geom_hline(yintercept = 0) +
    geom_line(data = nmf_result_ls$pure_rt[grepl("MS1", MSid), ]) +
    geom_line(data = nmf_result_ls$pure_rt[grepl("MS2", MSid), ], aes(y = -value_sc)) +
    facet_grid(rank ~ ., scales = "free_y") +
    theme_bw()

  ## Plot full iteration (mixed vs unmixed)



}
# Check feature 3008, 7797
    # Generating MS2 peaks list
    # Extracting ions signals
    # Generating 717 XICs
    # Build mixed matrix
    # Unmixing MS1Error in U[, i] : subscript out of bounds
{

  ms_xics_i[xic_label %in% rownames(ms_mixed_i$mixedmat),] %>% {
    ggplot(., aes(
      rtime,
      intensity,
      color = isolationWindowTargetMz,
      group = xic_label)
    ) +
      geom_line() +
      theme_bw()
  }

  ms_xics_i[, .N, by = .(xic_label)][order(-N),]
  ms_xics_i[xic_label %in% c("CP04497-1"), ]


  nmf_result_ls$pure_rt %>% {
    ggplot(., aes(
      rtime,
      value,
      color = MSid)
    ) +
      facet_grid(rank ~ ., scales = "free_y") +
      geom_line() +
      theme_bw()
  }

  mixed_mat_out %>% {
    ggplot(., aes(
      rtime,
      x,
      color = mslevel)
    ) +
      facet_grid(xic_label ~ .) +
      geom_line() +
      theme_bw()
  }
}