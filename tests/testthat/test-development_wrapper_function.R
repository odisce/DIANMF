testthat::skip("development script; test the wrapper function")

devtools::document()
devtools::load_all()
input_dir <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
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

xcms_obj <- DIANMF::detect_xcms_peaks(sequence_table = mzml_dt, params = params_ls)

features <- DIANMF::DIANMF.f(msexp = xcms_obj, d.out = FALSE,
                             sample_idx = 1,
                             MS2_ISOEACHL = T,
                             MS1MS2_L = F,
                             rank = 10,
                             maximumIteration = 200,
                             maxFBIteration = 100,
                             toleranceFB = 1e-05,
                             initialization_method = "nndsvd",
                             errors_print = FALSE,
                             method = "svds",
                             scan_rt_ext = 10,
                             min_distance = 5 )
