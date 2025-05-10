testthat::skip("development script; test the wrapper function")

devtools::document()
devtools::load_all()

## automatically detect path prefix
if(.Platform$OS.type == "unix") {
  pref <- "/spi"
  sysos <- "unix"
} else {
  pref <- "//fouet/spi"
  sysos <- "windows"
}
input_dir <- file.path(pref, "scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/")
input_targ <- file.path(pref, "scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Compound_data_corrected_rt.csv")
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
    binSize = 0.01,
    ppm = 7,
    maxFeatures = 500
  ),
  "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
)

cache_dir <- file.path(pref, "scidospi/06_Data/deconv_nmf", sysos)
cache_path <- file.path(cache_dir, "temp.rds")
profile_path <- file.path(pref, "scidospi/06_Data/deconv_nmf", sysos)
dir.create(dirname(cache_path))
if (file.exists(cache_path)) {
  xcms_obj <- readRDS(cache_path)
} else {
  xcms_obj <- DIANMF::detect_xcms_peaks(sequence_table = mzml_dt, params = params_ls)
  saveRDS(xcms_obj, cache_path)
}

{
  require(profvis)
  devtools::load_all()
  profile_save_path <- file.path(profile_path, format(Sys.time(), "%y%m%d-%H%M-profvis.html"))
  proffres <- profvis::profvis(
    expr = {
      features <- DIANMF.f(
        msexp = xcms_obj,
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
        scan_rt_ext = 10,
        min_distance = 4,
        featuresn = 2,
        nscans = 6,
        rt_method = "constant",
        clean_sources = TRUE,
        combineSpectra_arg = list(
          peaks = "intersect",
          ppm = 4,
          tolerance = 0.005,
          minProp = 0.05
        ),
        verbose = TRUE
      )
    }
  )
  htmlwidgets::saveWidget(proffres, profile_save_path)
}

cache_path <- file.path(cache_dir, "features.rds")
if (file.exists(cache_path)) {
  features <- readRDS(cache_path)
} else {
  data.table::setDTthreads(1)
  A <- Sys.time()
  # rt_range <- targ_dt[grepl("Scopolamine", Compound), (rt_sec + c(-60, +60))]
  # msexp <- xcms::filterRt(xcms_obj, rt_range)
  msexp <- xcms_obj
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
    scan_rt_ext = 10,
    min_distance = 4,
    featuresn = 2,
    nscans = 6,
    rt_method = "constant",
    clean_sources = TRUE,
    combineSpectra_arg = list(peaks = "intersect", ppm = 4, tolerance = 0.005, minProp = 0.05),
    verbose = TRUE
  )
  B <- Sys.time()
  difftime(B, A)
  # 2.28 hours for 3 files
  saveRDS(features, cache_path)
}

{
  ## Features summary
  ### Annotate features
  temp_ft <- get_feature_summary(features.l = features, max_method = "max_value")
  ft_match <- search_features(
    feature_dt = temp_ft,
    dt = targ_dt[, .(Compound, mz = mz_pos, rt = rt_sec)],
    rttol = 30,
    ppm = 10
  )
  temp_ft <- temp_ft[order(-into), ]

  ## Plot one sample diagnostic plot for one iteration
  save_path <- file.path(cache_dir, format(Sys.time(), "%y%m%d-%H%M"))
  dir.create(save_path)
  list.files(save_path, full.names = TRUE) %>% sapply(., unlink)
  for (feati in ft_match[, unique(featureid)]) {
    # feati <- "FT1501"
    timeA <- Sys.time()
    message(sprintf("Extracting: %s", feati), appendLF = FALSE)
    plot_out <- lapply(temp_ft[, unique(sample)], function(spli) {
      plot_feature(
        features.l = features,
        summary_dt = temp_ft,
        feature_id = feati,
        sample_index = spli,
        log2L = FALSE,
        max_method = "contribution",
        method = "best"
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

  feati <- ft_match[, unique(featureid)][1]
  ## Get best source coordinate for feati in all samples
  feat_coord <- get_feature_coord(
    features.l = features,
    feature_id = feati,
    sample_index = NULL,
    max_method = "contribution"
  )
  feat_coord[, .(featureid, sample, source, iteration, IsoWin)]

  ## Get best spectra for feati
  get_spectra(
    features.l = features,
    feature_id = feati,
    sample_index = 2,
    type = "pure",
    method = "best",
    max_method = "contribution"
  )
  ## get best pure elution profile
  get_elutionprofile(
    features.l = features,
    feature_id = feati,
    sample_index = 2,
    type = "pure",
    method = "best",
    max_method = "contribution"
  )
  ## Get best mixed spectra for feati
  get_spectra(
    features.l = features,
    feature_id = feati,
    sample_index = 2,
    type = "mixed",
    method = "best",
    max_method = "contribution"
  )
  ## get best mixed elution profile
  get_elutionprofile(
    features.l = features,
    feature_id = feati,
    sample_index = 2,
    type = "mixed",
    method = "best",
    max_method = "contribution"
  )
  ## get all sources for best iteration
  get_spectra(
    features.l = features,
    feature_id = feati,
    sample_index = 2,
    type = "pure",
    method = "all",
    max_method = "contribution"
  )
}