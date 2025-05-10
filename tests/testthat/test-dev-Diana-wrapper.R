testthat::skip("development script; test the wrapper function")

devtools::document()
devtools::load_all()

# remotes::install_git("http://132.166.53.123:8081/tools/r-packages/dianmf.git", ref = "dianmf_version2")
# this script is to process all Barbier mzML files
# library(DIANMF)
# input_dir <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
input_dir <- "~/mzml/"
mzml_dt <- DIANMF::prepare_mzMLfiles(input_dir = input_dir)
# mzml_dt <- mzml_dt[ class == 10,  ]
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
register(SnowParam(workers = 8))
# registered()
cache_path <- "~/dianmf_results/temp/xcms_obj.rds"
if (file.exists(cache_path)) {
  xcms_obj <- readRDS(cache_path)  
} else {
  xcms_obj <- DIANMF::detect_xcms_peaks(sequence_table = mzml_dt, params = params_ls)
  saveRDS(xcms_obj, cache_path)
}

# data.table::setDTthreads(1)
require(DIANMF)

features <- DIANMF.f(
  msexp = xcms_obj,
  dir_out = FALSE,
  sample_idx = NULL,
  MS2_ISOEACHL = TRUE,
  MS1MS2_L = FALSE,
  rank = 15,
  min_contrib = 0.6,
  maximumIteration = 200,
  maxFBIteration = 100,
  toleranceFB = 1e-05,
  initialization_method = "nndsvd",
  errors_print = FALSE,
  method = "svds",
  scan_rt_ext = 10,
  min_distance = 4,
  featuresn = NULL,
  nscans = 6,
  rt_method = "constant",
  clean_sources = TRUE,
  combineSpectra_arg = list(
    peaks = "intersect",
    ppm = 4,
    tolerance = 0.005,
    minProp = 0.05
  ),
  verbose = T
)

saveRDS(features, '~/dianmf_results/Results/features.rds')
