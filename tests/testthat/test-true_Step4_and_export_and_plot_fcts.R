test_that("step 4, export and plot functions tests", {

  input_dir <- system.file("extdata", package = "DIANMF")
  mzml_dt <- prepare_mzMLfiles(input_dir)
  mzml_seq <- create_seq(mzml_dt)

  params_ls <- list(
    "CentWaveParam" = xcms::CentWaveParam(
      ppm = 5,
      peakwidth = c(3, 15),
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
  xcms_obj <- detect_xcms_peaks(sequence_table = mzml_seq, params = params_ls)
  BiocParallel::register(SnowParam(workers = 1))
  
  features_sub <- DIANMF.f(
    msexp = xcms_obj,
    dir_out = FALSE,
    sample_idx = 1,
    MS2_ISOEACHL = TRUE,
    MS1MS2_L = TRUE,
    rank = 30,
    min_contrib = 0.6,
    maximumIteration = 200,
    maxFBIteration = 100,
    toleranceFB = 1e-05,
    initialization_method = "nndsvd",
    errors_print = FALSE,
    method = "svds",
    sparsityA = TRUE,
    scan_rt_ext = 10,
    min_distance = 4,
    featuresn = 2,
    nscans = 3,
    rt_method = "constant",
    clean_sources = TRUE,
    combineSpectra_arg = list(
      peaks = "union",
      ppm = 5,
      tolerance = 0.005
    ),
    verbose = T
  )
  
  expect_true(!is.null(features_sub))
  expect_true( all(c("PureFeatures", "ms1_features_peaks", "ms1_features") %in% names(features_sub[[1]])) )
  
  
  features <- DIANMF.f(
    msexp = xcms_obj,
    dir_out = FALSE,
    sample_idx = NULL,
    MS2_ISOEACHL = TRUE,
    MS1MS2_L = TRUE,
    rank = 30,
    min_contrib = 0.6,
    maximumIteration = 200,
    maxFBIteration = 100,
    toleranceFB = 1e-05,
    initialization_method = "nndsvd",
    errors_print = FALSE,
    method = "svds",
    sparsityA = TRUE,
    scan_rt_ext = 10,
    min_distance = 4,
    featuresn = 2,
    nscans = 3,
    rt_method = "constant",
    clean_sources = TRUE,
    combineSpectra_arg = list(
      peaks = "union",
      ppm = 5,
      tolerance = 0.005
      # minProp = 0.05
    ),
    verbose = T
  )
  
  temp_ft_max_value <- get_feature_summary(features.l = features, max_method = "max_value")
  temp_ft_contribution <- get_feature_summary(features.l = features, max_method = "contribution")
  
  expect_true( "data.table" %in% class(temp_ft_max_value) & "data.table" %in% class(temp_ft_contribution) )
  expect_true( nrow(temp_ft_max_value) > 0 & nrow(temp_ft_contribution) > 0 )
  
  feature_info <- get_feature_coord(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    max_method = "max_value" )
  
  expect_true( nrow(feature_info) == 1 )
  
  feature_elution_profiles_1 <- get_elutionprofile(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    type = c("pure", "mixed")[1],
    method = c("all", "best")[1],
    max_method = "max_value"
  )
  
  feature_elution_profiles_2 <- get_elutionprofile(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    type = c("pure", "mixed")[2],
    method = c("all", "best")[2],
    max_method = "contribution"
  )
  
  expect_true( nrow(feature_elution_profiles_1) > 1 & nrow(feature_elution_profiles_2) > 1 )
  expect_true( all(c("xic_label", "rtime", "value", "mslevel" ) %in% colnames(feature_elution_profiles_1)) )
  expect_true( all(c("xic_label", "rtime", "value", "mslevel" ) %in% colnames(feature_elution_profiles_2)) )
  
  
  feature_spectra_1 <- get_spectra(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    type = c("pure", "mixed")[1],
    method = c("all", "best")[1],
    max_method = "max_value"
  )
  
  feature_spectra_2 <- get_spectra(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    type = c("pure", "mixed")[2],
    method = c("all", "best")[2],
    max_method = "contribution"
  )
  
  expect_true( nrow(feature_spectra_1) > 0 & nrow(feature_spectra_2) > 0 )
  expect_true( all(c("xic_label", "IsoWin", "mslevel", "mz", "value", "rank" ) %in% colnames(feature_spectra_1)) )
  expect_true( all(c("xic_label", "IsoWin", "mslevel", "mz", "value", "rank" ) %in% colnames(feature_spectra_2)) )
  expect_true( nrow(feature_spectra_1)  >= nrow(feature_spectra_2) )
  
  
  chrom1_p <- plot_EluProfile(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    log2L = FALSE,
    type = c("pure", "mixed")[1],
    method = c("all", "best")[1],
    max_method
  )
  
  chrom2_p <- plot_EluProfile(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    log2L = TRUE,
    type = c("pure", "mixed")[1],
    method = c("all", "best")[1],
    max_method
  )
  
  expect_true( "ggplot" %in% class(chrom1_p) )
  expect_true( "ggplot" %in% class(chrom2_p) )
  
  
  spect1_p <- plot_Spectra(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    log2L = F,
    type = c("pure", "mixed")[1],
    method = c("all", "best")[1],
    max_method
  )
  
  spect2_p <- plot_Spectra(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    log2L = T,
    type = c("pure", "mixed")[1],
    method = c("all", "best")[1],
    max_method
  ) 
  
  expect_true( "ggplot" %in% class(spect1_p) )
  expect_true( "ggplot" %in% class(spect2_p) )
  
  
  feat_p <- plot_feature(
    features.l = features,
    summary_dt = temp_ft_max_value,
    feature_id = temp_ft_max_value$featureid[1],
    sample_index = 1,
    log2L = FALSE,
    max_method = "max_value",
    method = c("all", "best")[2]
  ) 
  
  expect_true( "ggplot" %in% class(feat_p) )
  
  peaks_p <- plot_xcms_peaks_range(summary_dt = temp_ft_max_value,
                        sample_index = 1,
                        iteration_index = 1)
    
  expect_true( "ggplot" %in% class(peaks_p) )
  
  rt_wind <- plot_rtWind_info(features.l = features, sample_index = 1, iteration_idx = 1, targets = NULL)
  rt_wind_p <- rt_wind$p_final 
  rt_wind_spect <- rt_wind$pure_spectra
  
  expect_true( "ggplot" %in% class(rt_wind_p) )
  expect_true( "data.table" %in%  class(rt_wind_spect) )
  expect_true( nrow(rt_wind_spect) > 0 )
  expect_true( all(c("xic_label", "rank", "value", "MSid", "scan_norm", "apex_val", "contribution", "mz", "rt", "mslevel") %in% colnames(rt_wind_spect)) )
  
  
  feat_spect_obj1 <- export_featureSpect(features.l = features,
                                        feature_id = temp_ft_max_value$featureid[1],
                                        sample_index = 1,
                                        type = c("pure", "mixed")[1],
                                        method = c("all", "best")[1],
                                        max_method = c("contribution", "max_value")[2] )
  
  expect_true( "Spectra" %in% class(feat_spect_obj1) )
  expect_true( length(mz(feat_spect_obj1)) == length(intensity(feat_spect_obj1))  )
  expect_true( length(mz(feat_spect_obj1)[[1]]) == length(intensity(feat_spect_obj1)[[1]])  )
  
  
  
  msSpectra <- exportMSSpectra(features.l = features,
                               sample_index = 1,
                               type = c("pure", "mixed")[1],
                               method = c("all", "best")[2],
                               max_method = c("contribution", "max_value")[2] )
  
  
  
  # MsBackendMgf::export(msSpectra, backend = MsBackendMgf(), file = "output.mgf")
  expect_true( "Spectra" %in% class(msSpectra) )
  expect_true( all(c("msLevel", "isolationWindowTargetMz", "name", "id" ) %in% spectraVariables(msSpectra))  )
  
  
  filtered_spect <- msSpectra[spectraData(msSpectra)$name == temp_ft_max_value$featureid[1]]
  expect_true( "Spectra" %in% class(filtered_spect) )
  
  ms2_spect <- filterMsLevel(filtered_spect, 2L)
  expect_true( "Spectra" %in% class(ms2_spect) )
  
  # ms2_spect_targMZ <- filterValues(ms2_spect, "isolationWindowTargetMz", 500) # I shoudn't fix 500
  # expect_true( "Spectra" %in% class(ms2_spect_targMZ) )
  
})