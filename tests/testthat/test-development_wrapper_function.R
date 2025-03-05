testthat::skip("development script; test the wrapper function")

devtools::document()
devtools::load_all()
input_dir <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
params_ls <- list(  
  "CentWaveParam" = CentWaveParam(
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
  "MergeNeighboringPeaksParam" = MergeNeighboringPeaksParam(
    expandRt = 2,
    expandMz = 0.001,
    ppm = 5,
    minProp = 0.75
  ),
  "ObiwarpParam" = ObiwarpParam(
    binSize = 0.05
  ),
  "PeakDensityParam" = PeakDensityParam(
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

features_mist <- DIANMF_f(input_dir, d.out = '~/DIA_NMF_R_package/Results4', 
                          sample_idx = 1, 
                          temp_saveL = T, MS2_ISOEACHL = T, MS1MS2_L = F,
                          scan_rt_ext = 10, min_distance = 5,
                          params = params_ls, plots = T )

