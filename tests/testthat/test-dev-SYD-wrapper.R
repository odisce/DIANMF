testthat::skip("development script; test the wrapper function")

devtools::document()
devtools::load_all()
input_dir <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
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

if (F) {
      msexp = xcms_obj
      dir_out = FALSE
      sample_idx = 1
      MS2_ISOEACHL = T
      MS1MS2_L = F
      rank = 10
      maximumIteration = 200
      maxFBIteration = 100
      toleranceFB = 1e-05
      initialization_method = "nndsvd"
      errors_print = FALSE
      method = "svds"
      scan_rt_ext = 5
      nscans = 6
      min_distance = 3
      featuresn = NULL
      verbose = TRUE
      rt_method = "constant"
      combineSpectra_arg = list(peaks = "intersect", ppm = 4, tolerance = 0.005, minProp = 0.05)
}

devtools::load_all()
A <- Sys.time()

    features <- DIANMF.f(
      msexp = xcms_obj,
      dir_out = FALSE,
      sample_idx = 1,
      MS2_ISOEACHL = TRUE,
      MS1MS2_L = FALSE,
      rank = 20,
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
# Check feature 3008, 7797
    # Generating MS2 peaks list
    # Extracting ions signals
    # Generating 717 XICs
    # Build mixed matrix
    # Unmixing MS1Error in U[, i] : subscript out of bounds
{
  require(ggplot2)
  ggplot(ms1_peaks[sample == 1 & is_filled == 0,]) +
    geom_linerange(aes(y = mz, x = rt, xmin = rtmin, xmax = rtmax)) +
    geom_linerange(data = ms1_peaks_i[sample == 1 & is_filled == 0,], aes(y = mz, x = rt, xmin = rtmin, xmax = rtmax, color = peakfull)) + 
    theme_bw()

    ggplot(pure_rt_good, aes(rtime, value, color = rank, group = rank)) +
      geom_line() +
      facet_grid(rank~MSid) +
      theme_bw()
  
}