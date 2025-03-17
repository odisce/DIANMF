testthat::skip("development script; test the wrapper function")

{
  devtools::document()
  devtools::load_all()
  input_dir <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
  input_targ <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Compound_data_corrected_rt.csv"
  targ_dt <- fread(input_targ)
  mzml_dt <- prepare_mzMLfiles(input_dir = input_dir)
  mzml_dt <- mzml_dt[1:3, ]
  params_ls <- list(  
    "CentWaveParam" = xcms::CentWaveParam(
      ppm = 6,
      peakwidth = c(6, 30),
      snthresh = 0,
      prefilter = c(5, 4000),
      mzCenterFun = "wMeanApex3",
      integrate = 2,
      mzdiff = -0.001,
      noise = 2000,
      firstBaselineCheck = FALSE
    ),
    "MergeNeighboringPeaksParam" = xcms::MergeNeighboringPeaksParam(
      expandRt = 2,
      expandMz = 0.001,
      ppm = 5,
      minProp = 0.75
    ),
    "ObiwarpParam" = xcms::ObiwarpParam(
      binSize = 0.05
    ),
    "PeakDensityParam" = xcms::PeakDensityParam(
      sampleGroups = NA,
      bw = 15,
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
        featuresn = 3
      )
    }
  )
  htmlwidgets::saveWidget(proffres, profile_save_path)
}

devtools::load_all()
data.table::setDTthreads(1)
A <- Sys.time()
    rt_range <- targ_dt[grepl("Scopolamine", Compound), (rt_sec + c(-60, +60))]
    features <- DIANMF.f(
      msexp = xcms::filterRt(xcms_obj, rt_range),
      dir_out = FALSE,
      sample_idx = 1,
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
      scan_rt_ext = 30,
      min_distance = 4,
      featuresn = NULL,
      rt_method = "constant",
      combineSpectra_arg = list(peaks = "intersect", ppm = 4, tolerance = 0.005, minProp = 0.05),
      verbose = TRUE
    )

B <- Sys.time()

difftime(B, A)

{
  ## Get pure rt and spectra for feature n
  features[[1]]$ms1_features
  featuresel <- "FT1834"
  lapply(features, function(x) {
    x$ms1_features[featureid == featuresel,] %>% na.omit()
  })

  feature_dt <- data.table()
  for (spli in seq_len(length(features))) {
    for (iteri in seq_len(length(features[[spli]]$PureFeatures))) {
      # spli <- 1 ; iteri <- 1
      if (is.null(features[[spli]]$PureFeatures[[iteri]])) {
        next
      }
      # print(sprintf("spl: %i, feat: %i", spli, iteri))
      # names(features[[spli]]$PureFeatures[[iteri]])
      out_dt <- features[[spli]]$PureFeatures[[iteri]]$ms_info[!is.na(source), .(featureid, source, sample = spli, iteration = iteri)]
      # print(nrow(out_dt))
      feature_dt <- rbind(feature_dt, out_dt)
      rm(out_dt)
    }
  }
  

  ## Features summary

  ## Plot one sample diagnostic plot for one iteration
  names(features[[1]])
  names(features[[1]]$PureFeatures[[1]])
  features[[1]]$PureFeatures[[3]]$ms2_info[!is.na(source), ]

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

  
}