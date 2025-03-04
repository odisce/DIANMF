testthat::skip("development script to implement xcms sequence objects")

## Improve XCMS parameters
## Run on Barbier data

library(MSnbase)
library(data.table)
library(BiocParallel)
library(magrittr)
library(ggplot2)
library(ggpubr)
devtools::load_all()
input_dir <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
input_target <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Compound_data_ids.csv"
input_targ <- fread(input_target)

## Create sequence
input_files <- list.files(input_dir, pattern = ".mzml", full.names = TRUE, ignore.case = TRUE)
mzml_dt <- data.table(
  "mzml_path" = input_files,
  "file_name" = tools::file_path_sans_ext(basename(input_files)),
  "creation_date" = file.info(input_files)$ctime %>% as.POSIXct(),
  "class" = gsub("^.*-([0-9]{1,3}(,[0-9]{1,3})?)ng.*$", "\\1", basename(input_files)) %>%
    gsub(",", ".", .) %>%
    as.numeric()
)
mzml_dt[order(creation_date), InjectionOrder := seq_len(.N)]

# xcms parameters 
params_ls <- list(  
  "CentWaveParam" = CentWaveParam(
    ppm = 6,
    peakwidth = c(6, 30),
    snthresh = 0,
    prefilter = c(5, 4000),
    mzCenterFun = "wMeanApex3",
    integrate = 2,
    mzdiff = -0.0003,
    noise = 2000,
    firstBaselineCheck = FALSE
  ),
  "MergeNeighboringPeaksParam" = MergeNeighboringPeaksParam(
    expandRt = 2,
    expandMz = 0,
    ppm = 1,
    minProp = 0.75
  ),
  "ObiwarpParam" = ObiwarpParam(
    binSize = 0.05
  ),
  "PeakDensityParam" = PeakDensityParam(
    sampleGroups = NA,
    bw = 15,
    minFraction = 0,
    minSamples = 1,
    binSize = 0.005,
    ppm = 5,
    maxFeatures = 500
  ),
  "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
)


# register(SerialParam())
register(MulticoreParam(3))
bpparam()

## Detect LCpeaks ----------------------------
temp_saveL <- T
subset_rt <- 297 + c(-10, +10) 
known_compounds <- input_targ[rt_sec %between% subset_rt ]

target_mz <- 304.1543 
# xcms_obj <- detect_xcms_peaks(mzml_dt[, lapply(.SD, head, 1), by = .(class)], params_ls)
if (temp_saveL) {
  save_path <- "./temp/data/"
  dir.create(save_path, recursive = TRUE)
  xcms_obj_path <- file.path(save_path, "xcms_obj.rds")
  if (file.exists(xcms_obj_path)) {
    xcms_obj <- readRDS(xcms_obj_path)
  } else {
    xcms_obj <- detect_xcms_peaks(
      sequence_table = mzml_dt[, lapply(.SD, head, 1), by = .(class)],
      params = params_ls,
      rt_range = subset_rt
    )
    saveRDS(xcms_obj, file = xcms_obj_path)
  }
}

## Order peak targets by intensities ----------------------------
ms1_peaks <- extract_xcms_peaks(xcms_obj, 1, TRUE)

## Iterate over peaks for one sample ----------------------------
ms1_peaks[, cpnt := as.integer(NA)]
cpnt <- 1
# for (i in seq_len(ms1_peaks[, .N])) {
# i <- 1
i <- search_target_peak <- ms1_peaks[, which(mz %between% (target_mz + c(-0.001, +0.001)))]

res <- lapply(1:nrow(known_compounds), function(h){
  # h = 2
  i_known_compounds <-  ms1_peaks[, which(mz %between% (known_compounds$mz_pos[[h]] + c(-0.001, +0.001)))]
  peak_i <- ms1_peaks[i_known_compounds, ]
  peak_i$compound <- known_compounds$Compound[[h]]
  peak_i
})
names(res) <- known_compounds$Compound
res <- do.call(rbind, res)

useapexL <- T
peak_i <- ms1_peaks[i, ]
if(useapexL) {
  ## use apex window
  scan_rt <- 10
  rt_range <- peak_i[, (rt + c(-scan_rt-1, +scan_rt+1))]
  peaks_i <- ms1_peaks[rtmin <= rt_range[2] & rtmax >= rt_range[1], ]
} else {
  ## full peak range
  rt_range <- peak_i[, c(rtmin, rtmax)]
  peaks_i <- ms1_peaks[rtmin <= peak_i[, rtmax] & rtmax >= peak_i[, rtmin], ]
}

## flag peaks partially, fully or apex inside the window
peaks_i[, peakfull := ifelse(
  (rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
  ifelse(rt %between% rt_range, "apex", "partial")
)]
{ ## plot peaks in range
  ggplot(peaks_i, aes(rt, mz, group = peakid)) +
    geom_pointrange(aes(x = rt, xmin = rtmin, xmax = rtmax, color = peakfull)) +
    geom_pointrange(data = res, aes(x = rt, xmin = rtmin, xmax = rtmax), color = "black") +
    geom_text(data = res, 
              aes(x = rt, y = mz, label = compound),
              color = 'black', vjust = -0.5, size = 3.5) +
    # geom_pointrange(data = peak_i, aes(x = rt, xmin = rtmin, xmax = rtmax), color = "gray") +
    geom_vline(xintercept = rt_range, linetype = 2, color = "red") +
    theme_bw() +
    labs(
      title = "MS1 peak in range",
      caption = sprintf(
        "
            Total MS1 peaks: %i\n
            Peak include in range: %i
          ",
        peaks_i[, .N],
        peaks_i[peakfull == TRUE, .N]
      )
    )
}

## extract data
# raw_dt <- xcms::filterRt(xcms_obj, rt_range) %>%
#   xcms::filterFile(., 1) %>%
#   # xcms::filterMsLevel(., 1) %>%
#   get_spectra_values()
# ## Normalize time
# time_dic <- raw_dt[, .(rtime, msLevel, isolationWindowTargetMz)] %>% unique()
# time_dic[order(rtime), scan_norm := seq_len(.N), by = .(msLevel, isolationWindowTargetMz)]
# raw_dt <- merge(
#   raw_dt,
#   time_dic,
#   by = c("rtime", "msLevel", "isolationWindowTargetMz")
# )
# alignement
raw_dt <- xcms::filterRt(xcms_obj, rt_range) %>%
  xcms::filterFile(., 1) %>%
  # xcms::filterMsLevel(., 1) %>%
  get_spectra_values()
## Normalize time
time_dic <- raw_dt[, .(rtime, msLevel, isolationWindowTargetMz)] %>% unique()
time_dic_ms1 <- time_dic[order(rtime), ][msLevel == 1, ]
time_dic_ms1[, scan_norm := seq_len(.N)]
ms2_rtime <- time_dic_ms1[, {
  rt_ref <- rtime
  time_dic[order(abs(rtime - rt_ref)), ][, lapply(.SD, head, 1), by = .(isolationWindowTargetMz)]
}, by = .(scan_norm)]
time_dic <- ms2_rtime[scan_norm %between% (range(scan_norm) + c(+1, -1)), ]
time_dic[order(rtime), scan_norm := seq_len(.N), by = .(msLevel, isolationWindowTargetMz)]
dcast(time_dic, scan_norm ~ isolationWindowTargetMz, value.var = 'rtime')
raw_dt[, .N, by = .(rtime, msLevel, isolationWindowTargetMz)]
raw_dt <- merge(
  raw_dt,
  time_dic,
  by = c("rtime", "msLevel", "isolationWindowTargetMz"),
  allow.cartesian = T
)


## build MS1 xics from peak list
xic_dt_ms1 <- peaks_i[, {
  mzrange <- c(mzmin, mzmax)
  raw_dt[msLevel == 1 & mz %between% mzrange, .(mz, scan_norm, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy)]
}, by = .(peakid, peakfull)]
xic_dt_ms1[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
xic_dt_ms1 <- xic_dt_ms1[xic_label %in% xic_dt_ms1[, .N, by = xic_label][N >= 4, xic_label]]  # filter MS1 eics

# ## build MS1 xics from raw data
# ### Combine all MS1 spectra
# MS1_peaklist <- xcms::filterRt(xcms_obj, rt_range) %>%
#   xcms::filterFile(., 1) %>%
#   xcms::filterMsLevel(., 1) %>%
#   xcms::spectra() %>%
#   Spectra::combineSpectra(ppm = 3, tolerance = 0.001, minProp = 0.001) %>%
#   Spectra::asDataFrame() %>%
#   data.table::as.data.table() %>%
#   {.[, .(msLevel, mz, rtime, intensity)]}
# MS1_peaklist[, xic_label := paste0(seq_len(.N), "-", 1)]
# ### Extract XICs
# xic_dt_ms1 <- MS1_peaklist[, {
#   mzrange <- PpmRange(mz, 7) 
#   raw_dt[msLevel == 1 & mz %between% mzrange, .(mz, rtime, intensity, msLevel)]
# }, by = .(xic_label)]
# xic_dt_ms1 <- xic_dt_ms1[xic_label %in% xic_dt_ms1[, .N, by = xic_label][N >= 4, xic_label]]
# xic_dt_ms1 <- merge(
#   xic_dt_ms1,
#   time_dic,
#   by = c("rtime", "msLevel")
# )

## build MS2 xics from raw data
### Combine all MS2 spectra

MS2_ISOEACHL <- T  # from each isolation window separately
if (MS2_ISOEACHL) {
  isowin <- xcms_obj %>%
    xcms::filterFile(., 1) %>%
    xcms::filterRt(., rt_range) %>%
    xcms::filterMsLevel(., 2L) %>%
    xcms::spectra() %>%
    xcms::isolationWindowTargetMz() %>%
    as.numeric() %>% unique()
  MS2_peaklist_l <- list()
  for (n in seq_len(length(isowin))) {
    i <- isowin[[n]]
    MS2_peaklist_l[[n]] <- xcms::filterRt(xcms_obj, rt_range + c(+1, -1)) %>%
      xcms::filterFile(., 1) %>%
      xcms::filterMsLevel(., 2) %>%
      xcms::spectra() %>%
      Spectra::filterIsolationWindow(i) %>%
      Spectra::combineSpectra(ppm = 3, tolerance = 0.005, minProp = 0.001) %>%
      Spectra::asDataFrame() %>%
      data.table::as.data.table() %>%
      {.[, .(msLevel, mz, rtime, intensity, isolationWindowTargetMz)]}
  }
  MS2_peaklist <- data.table::rbindlist(MS2_peaklist_l, idcol = "IsoWin")
  MS2_peaklist[, xic_label := paste0(seq_len(.N), "-", IsoWin, "-", 2)]
} else {
  MS2_peaklist <- xcms::filterRt(xcms_obj, rt_range) %>%
    xcms::filterFile(., 1) %>%
    xcms::filterMsLevel(., 2) %>%
    xcms::spectra() %>%
    Spectra::combineSpectra(ppm = 5, tolerance = 0.008, minProp = 0.001) %>%  # (ppm = 5, tolerance = 0.008, minProp = 0.01)
    Spectra::asDataFrame() %>%
    data.table::as.data.table() %>%
    {.[, .(msLevel, mz, rtime, intensity)]}
  MS2_peaklist[, xic_label := paste0(seq_len(.N), "-", 2)]
}

### Extract XICs
xic_dt_ms2 <- MS2_peaklist[, {
  mzrange <- PpmRange(mz, 7)
  raw_dt[msLevel == 2 & mz %between% mzrange, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz,
                                                 lowMZ, highMZ, collisionEnergy)]
}, by = .(xic_label)]
xic_dt_ms2 <- xic_dt_ms2[xic_label %in% xic_dt_ms2[, .N, by = xic_label][N >= 4, xic_label]]
xic_dt_ms2 <- merge(
  xic_dt_ms2,
  time_dic,
  by = c("rtime", "msLevel", "isolationWindowTargetMz"),
  allow.cartesian = T
)
xic_dt_ms2[, xic_label := paste0(xic_label, "-", isolationWindowTargetMz)]


# {
#   ## plot MS2
#   {
#     ggplot(
#       xic_dt_ms2[isolationWindowTargetMz %in% unique(isolationWindowTargetMz)[1:2],],
#       aes(mz, intensity, color = as.factor(isolationWindowTargetMz))
#     ) +
#       geom_linerange(aes(ymin = 0, ymax = intensity, text = xic_label)) +
#       theme_bw()
#   } %>%
#     plotly::ggplotly()
# }
{
  # plot xics
  plot_peak_map <- ggplot(xic_dt_ms1, aes(rtime, mz, group = peakid)) +
    geom_line(aes(color = peakfull)) +
    geom_line(data = xic_dt_ms1[peakid == peak_i[, peakid], ], color = "black") +
    facet_grid(msLevel ~ .) +
    theme_bw()

  xic_dt_norm <- xic_dt_ms1%>%
    group_by(peakid) %>%
    mutate(intensity_norm = intensity/sum(intensity))
  xic_dt_norm <- as.data.table(xic_dt_norm)
  xic_dt_norm_sub <- xic_dt_norm[peakid %in% res[, peakid], ]
  xic_dt_norm_sub <- merge(xic_dt_norm_sub, res, by = 'peakid')
  xic_dt_norm_sub <- xic_dt_norm_sub %>%
    group_by(peakid) %>%
    arrange(-intensity_norm) %>%
    slice(1)
  
  plot_xics <- ggplot(xic_dt_norm, aes(rtime, intensity, group = peakid)) +
    geom_line(aes(color = peakfull)) +
    geom_line(data = xic_dt_norm[peakid %in% res[, peakid], ],  aes( rtime, intensity), color = "black", linewidth = 0.8) +
    geom_text(data = xic_dt_norm_sub,
              aes(x = rt, y = intensity, label = compound),
              color = 'red', angle = 0, vjust = 0, size = 3) +
    geom_point() +
    facet_grid(msLevel ~ .) +
    theme_bw() +
    guides(color = "none") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
  
  # plot_xics <- ggplot(xic_dt_ms1, aes(rtime, intensity, group = xic_label, color = xic_label)) +
  #   geom_line() +
  #   geom_point() +
  #   theme_bw() +
  #   guides(color = "none") +
  #   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
  
  plot_xics_ms2 <- ggplot(xic_dt_ms2, aes(scan_norm, intensity, group = xic_label, color = "red") ) +
    geom_line(data = xic_dt_ms1, aes( y = intensity, color = "black")) +
    geom_line() + # aes(color = as.factor(isolationWindowTargetMz))
    theme_bw() +
    guides(color = "none") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

  plot_spectrum <- ggplot(peaks_i, aes(mz, into, group = peakid)) +
    geom_linerange(aes(color = peakfull, ymin = 0, ymax = into)) +
    geom_linerange(data = peak_i, color = "black", aes(color = peakfull, ymin = 0, ymax = into)) +
    facet_grid("1" ~ .) +
    theme_bw() +
    guides(color = "none") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
  
  ggpubr::ggarrange(
    plot_xics,
    plot_xics_ms2,
    nrow = 2,
    align = 'hv',
    heights = c(1,3)
  )
  # plot_peak_map
  # plot_spectrum
}


## build matrix (only on xcms peaks for now)
MS1MS2_L <- T
# ms1 matrix from xcms eics
ms1_mixeddt <- dcast(xic_dt_ms1[msLevel == 1, ], xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
ms1_mixedmat <- as.matrix(ms1_mixeddt, rownames = TRUE)
row_filter_ms1 <- apply(ms1_mixedmat, 1, has_four_consecutive_non_zero)
ms1_mixedmat <- ms1_mixedmat[row_filter_ms1, , drop = FALSE]

# deleted_ms1_part <- as.matrix(ms1_mixeddt, rownames = TRUE)
# deleted_ms1_part <- ms1_mixeddt[!row_filter_ms1, , drop = FALSE]

ms1_infos <- xic_dt_ms1[msLevel == 1, ][, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]
ms1_infos$isolationWindowLowerMz <- NA
ms1_infos$isolationWindowUpperMz <- NA
ms1_infos$lowMZ <- NA
ms1_infos$highMZ <- NA

# # ms1 matrix from raw data
# ms1_mixeddt <- dcast(xic_dt_ms1, xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
# ms1_mixedmat <- as.matrix(ms1_mixeddt, rownames = TRUE)
# ms1_infos <- xic_dt_ms1[msLevel == 1, ][, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]

ms2_mixeddt <- dcast(xic_dt_ms2, xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
ms2_mixedmat <- as.matrix(ms2_mixeddt, rownames = TRUE)
row_filter_ms2 <- apply(ms2_mixedmat, 1, has_four_consecutive_non_zero)
ms2_mixedmat <- ms2_mixedmat[row_filter_ms2, , drop = FALSE]
ms2_infos <- xic_dt_ms2[, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz, lowMZ, highMZ)]

if(MS1MS2_L) {
  mixedmat <- rbind(
    ms1_mixedmat,
    ms2_mixedmat )
  
  ms_infos <- rbind(
    ms1_infos,
    ms2_infos
  )
  dim(mixedmat)
} else {
  mixedmat <- ms1_mixedmat
  ms_infos <- ms1_infos
}
# REMOVE_COLINEARITY_L <- FALSE

## NMF
rank <- 15

# ms_mixed <- mixedmat
# ms_mixed <- reshape2::melt(ms_mixed)
# colnames(ms_mixed) <- c('mz_value', 'rt', 'intensity')
# ggplot(ms_mixed, aes(rt, intensity, color = as.factor(mz_value) )) +
#   geom_line() +
#   geom_point() +
#   theme_bw() +
#   guides(color = "none") +
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

#-----------------------------------------------------------------------------

ngmcas_res <- nGMCAs(
  X.m = mixedmat,
  rank = rank,
  maximumIteration = 200,
  maxFBIteration = 100,
  toleranceFB = 1e-05,
  initialization_method = "nndsvd",
  errors_print = FALSE,
  method = "svsd"
)
{
  
  if(MS1MS2_L == F){ # just MS1
    
    pure_rt_ms1 <- melt(ngmcas_res$A) %>% as.data.table()
    setnames(pure_rt_ms1, c("rank", "scan_norm", "value"))
    
    S <- ngmcas_res$S
    pure_mz_ms1 <- melt(S) %>% as.data.table()
    setnames(pure_mz_ms1, c("xic_label", "rank", "value"))
    pure_mz_ms1 <- merge(ms1_infos, pure_mz_ms1, by = "xic_label")
    
    ggpubr::ggarrange(
      ggpubr::ggarrange(
        plot_xics,
        plot_spectrum,
        align = "hv",
        ncol = 2
      ),
      ggpubr::ggarrange(
        ggplot(pure_rt_ms1, aes(scan_norm, value, group = as.factor(rank))) +
          geom_line() +
          facet_grid(rank ~ ., scales = "free_y") +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
        ggplot() +
          geom_hline(yintercept = 0) +
          geom_linerange(data = pure_mz_ms1[msLevel == 1, ], aes(mz, value, ymin = 0, ymax = value)) +
          # geom_linerange(data = pure_mz_ms2[msLevel == 2, ], aes(mz, -value, ymin = 0, ymax = -value), color = "red") +
          facet_grid(rank ~ .) +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
          guides(color = "none"),
        align = "hv",
        ncol = 2
      ),
      nrow = 2,
      heights = c(1, 2)
    )
  } else { # MS1 and MS2 combined
    
    # extract the MS1 data
    S_ms1 <- ngmcas_res$S[1:dim(ms1_mixedmat)[1], ]
    pure_rt_ms1 <- melt(ngmcas_res$A) %>% as.data.table()
    setnames(pure_rt_ms1, c("rank", "scan_norm", "value"))
    pure_mz_ms1 <- melt(S_ms1) %>% as.data.table()
    setnames(pure_mz_ms1, c("xic_label", "rank", "value"))
    pure_mz_ms1 <- merge(ms1_infos, pure_mz_ms1, by = "xic_label")
    
    # extract the MS2 data
    S_ms2 <- ngmcas_res$S[(dim(ms1_mixedmat)[1]+1):( (dim(ms1_mixedmat)[1]) + dim(ms2_mixedmat)[1]), ]
    pure_rt_ms2 <- melt(ngmcas_res$A) %>% as.data.table()
    setnames(pure_rt_ms2, c("rank", "scan_norm", "value"))
    pure_mz_ms2 <- melt(S_ms2) %>% as.data.table()
    setnames(pure_mz_ms2, c("xic_label", "rank", "value"))
    pure_mz_ms2 <- merge(ms2_infos, pure_mz_ms2, by = "xic_label", allow.cartesian=TRUE)
   
    ggpubr::ggarrange(
      ggpubr::ggarrange(
        plot_xics,
        plot_spectrum,
        align = "hv",
        ncol = 2
      ),
      ggpubr::ggarrange(
        ggplot(pure_rt_ms1, aes(scan_norm, value, group = as.factor(rank))) +
          geom_line() +
          facet_grid(rank ~ ., scales = "free_y") +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
        ggplot() +
          geom_hline(yintercept = 0) +
          geom_linerange(data = pure_mz_ms1[msLevel == 1, ], aes(mz, value, ymin = 0, ymax = value)) +
          geom_linerange(data = pure_mz_ms2[msLevel == 2, ], aes(mz, -value, ymin = 0, ymax = -value), color = "red") +
          facet_grid(rank ~ .) +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
          guides(color = "none"),
        align = "hv",
        ncol = 2
      ),
      nrow = 2,
      heights = c(1, 2)
    )
  }
  
}

#--------------------------------------------------------------------------------
# prepare MS1 pure sources and match them
{
  NewDb_ms1 <- readRDS("//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/REF_database/MS1_database.rds")
  #   
  # W <- ngmcas_res$S
  # rownames(W) <- pure_mz[rank == 1]$mz
  # H <- ngmcas_res$A
  # colnames(H) <- pure_rt[ rank == 1, ]$scan_norm
  # ms1_pure_sources <- pure_sources.f(W = W, H = H, ms_type = 'sum')
  
  subset_rt <- 297 + c(-10, +10)
  exist_comp <- input_targ[rt_sec %between% (subset_rt)]
  ms1_mz_tolerance <- 0.007

  # library(future.apply)
  # plan(multisession, workers = 4)

  res <- lapply(unique(pure_mz_ms1$rank), function(s){
    # s = 1
    print(s)
    ms1_spectrum <- pure_mz_ms1[rank == s & msLevel == 1]
    colnames(ms1_spectrum)[c(4,10)] <- c( "mz_value", "intensity")
    ms1_spectrum <- ms1_spectrum[ms1_spectrum$intensity > 0, ]
    
    if( nrow(ms1_spectrum) > 1 & !is.null(ms1_spectrum) ){
      ms1_spectrum$intensity <- ms1_spectrum$intensity / max(ms1_spectrum$intensity)
      ms1_spectrum$mz_value <- gsub("\\.\\d$", "",  ms1_spectrum$mz_value) # Remove trailing '.X' where X is any digit
      ms1_spectrum$mz_value <- as.numeric( ms1_spectrum$mz_value)
      ms1_spectrum <- as.data.frame(ms1_spectrum)
      
      # ms1_spectrum <- ms1_pure_sources[[s]]$source_spect
      # ms1_spectrum <- ms1_spectrum[ms1_spectrum$intensity > 0, ]
      # ms1_spectrum$intensity <- ms1_spectrum$intensity / max(ms1_spectrum$intensity)
      # ms1_spectrum <- as.data.frame(ms1_spectrum)
      # ms1_spectrum$mz_value <- gsub("\\.\\d$", "",  ms1_spectrum$mz_value) # Remove trailing '.X' where X is any digit
      # ms1_spectrum$mz_value <- as.numeric( ms1_spectrum$mz_value)
      
      res_sub1 <- lapply(unique(exist_comp$chemical_id), function(c){ #   names(NewDb_ms1)
        # c <- "334"
        comp <- exist_comp[ exist_comp$chemical_id == c, ]$Compound
        c <- as.character(c)
        print(c)
        ms1_lib_spect <- NewDb_ms1[[c]]
        ms1_lib_spect <- list(ms1_lib_spect)
        if( !is.null(ms1_lib_spect) ){
          ms1_score <- lapply(1:length(ms1_lib_spect), function(j){
            # j <- 1
            # print(j)
            ms1_lib_spectrum <- ms1_lib_spect[[j]]
            # ms1_lib_spectrum <- ms1_lib_spectrum[, c(4,3)]
            colnames(ms1_lib_spectrum)[c(6,7)] <- c('mz_value', 'intensity')
            ms1_lib_spectrum <- as.data.frame(ms1_lib_spectrum)
            # ms1_lib_spectrum$intensity <- ms1_lib_spectrum$intensity / max(ms1_lib_spectrum$intensity)
            score1 <- round(GetSimpleDotProductSimilarity(ms1_spectrum, ms1_lib_spectrum, bin = ms1_mz_tolerance), 3)
            score2 <-  round(getReverseSearchingSimilarity(measuredSpectra = ms1_spectrum, librarySpectra = ms1_lib_spectrum, bin = ms1_mz_tolerance), 3)
            score3 <- round(GetPresenceSimilarity(ms1_spectrum, ms1_lib_spectrum, bin = ms1_mz_tolerance), 3) 
            score <- ( 0.1*score1 + 0.4*score2 + 0.5*score3 )
            
            data.table(
              "source" = s,
              "name" = comp, 
              "chemical_id" = as.numeric(c),
              "lib_spect_idx" = j,
              "dot_prod" = score1,
              "rev_prod" = score2,
              "presence" = score3,
              "total_score" = score )
          })
          ms1_score <- do.call(rbind, ms1_score)
          ms1_score <- ms1_score[ms1_score$total_score > 0, ]
          ms1_score
          if( nrow(ms1_score) == 0 ){
            ms1_score <- NULL
          }else{
            lib_indx <- which.max(ms1_score$total_score)
            ms1_score <- ms1_score[lib_indx, ]
          }
        }
      })
      names(res_sub1) <- exist_comp$Compound
      res_sub1 <- do.call(rbind, res_sub1)
      
      d.plot <- file.path('~/DIA_NMF_R_package/Results2/ms1/')
      if (!dir.exists(d.plot)) { dir.create(d.plot, recursive = TRUE) }
      if(max(res_sub1$total_score) > 0){
        res_sub2 <- res_sub1 %>% 
          arrange(desc(total_score)) %>%
          slice_head(n = 1)
        
        plots_list <- list()
        for(i in 1:nrow(res_sub2)){
          # i <- 1
          r <- res_sub2[i, ]
          ms1_lib_spectrum <- NewDb_ms1[[ as.character(r$chemical_id)]]
          colnames(ms1_lib_spectrum)[c(6,7)] <- c('mz_value', 'intensity')
          ms1_lib_spectrum <- as.data.frame(ms1_lib_spectrum)

          p1 <- ggplot() +
            geom_linerange(data = ms1_spectrum, 
                           aes(x = mz_value, ymin = 0, ymax = intensity), color = 'blue') +
            geom_linerange(data = ms1_lib_spectrum, aes(x = mz_value, ymin = -intensity, ymax = 0), color = 'red') +
            guides(color = "none") + 
            labs(x = "mz",
                 y = "Intensity",
                 caption = paste('workflow =', "DIANMF", ':',
                                 'DP =',  r$dot_prod, 
                                 '| IDP =', r$rev_prod, 
                                 '| FP =', r$presence, 
                                 '| Name =', r$name)) +
            theme_bw(base_size = 14) +
            theme(
              plot.title = element_text( hjust = 0.5, size = 16),
              axis.title = element_text(),
              axis.text = element_text(size = 12),
              plot.caption = element_text(hjust = 0, face = "italic", size = 10),
              legend.position = "none"   )
          
          plots_list[[i]] <- p1
        }
        p_all <- patchwork::wrap_plots(plots_list, nrow = 1)
        
        eics_mat <- pure_rt_ms1[rank == s, ]
        p1 <- ggplot(eics_mat, aes(x = scan_norm, y = value)) +
          geom_line() +
          geom_point() +
          xlim(min(eics_mat$scan_norm), max(eics_mat$scan_norm)) +
          guides(color = 'none') +
          theme_minimal()
        p <- list(
          'spect' = p_all,
          'eics' = p1  )
        p <- patchwork::wrap_plots(p, nrow = 2)
      
        ggsave(paste0(d.plot, '/', 'source_', r$source, "_", r$name, '_ms1.png'), p, w=15, h = 15, dpi = 300)
      }
      
      return(res_sub1)
    }
  })
  ress <- res
  names(ress) <- paste0("source_", seq(1,rank))
  ress <- do.call(rbind, ress)
  
  ress2 <- ress %>% 
    group_by(source) %>%
    arrange(desc(total_score)) %>%
    slice_head(n = 1)
  saveRDS(ress2, "~/DIA_NMF_R_package/Results2/ms1/ms1_scores.rds")
}

# contribution
{
  pure_mz_ms1 <- pure_mz_ms1[ msLevel == 1, ]
  pure_mz_ms1 <- pure_mz_ms1 %>%
    group_by(rank) %>%
    mutate(total_value = sum(value, na.rm = TRUE),
           contribution = value / total_value) %>%
    ungroup()
  
  mz_tolerance <- 0.01  
  
  pure_mz_ms1 <- pure_mz_ms1 %>%
    group_by(rank) %>%
    mutate(value_norm = value/max(value)) %>%
    ungroup()
  
  pure_mz_ms1 <- as.data.table(pure_mz_ms1)
  
  pure_mz_ms1_sub <- lapply(1:nrow(known_compounds), function(s){
    mz_subset <- known_compounds[s, ]$mz_pos + c(-0.001, +0.001)
    ions <- pure_mz_ms1[ mz %between% mz_subset, ]
    ions <- ions[ value > 0,  ]
    ions$compound <- known_compounds[s, ]$Compound
    return(ions)
  })
  pure_mz_ms1_sub <- do.call(rbind, pure_mz_ms1_sub)
  
  p <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_linerange(data = pure_mz_ms1, aes(mz, value, ymin = 0, ymax = value)) +
    geom_linerange(data = pure_mz_ms1_sub, aes(mz, value, ymin = 0, ymax = value, color = "red") ) +
    geom_text(data = pure_mz_ms1_sub, 
              aes(x = mz, y = value, label = paste(compound, round(contribution,3) )),
              color = 'red', vjust = 0.1, size = 2.5) +
    facet_wrap(rank ~ ., ncol = 2) +
    theme_bw() +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    guides(color = "none")
  
  d.plot <- file.path('~/DIA_NMF_R_package/Results2/ms1/')
  ggsave(paste0(d.plot, '/', 'ms1_contribution.png'), p, w=10, h = 10, dpi = 300)
  
  p
}

# ------------------------------------------------------------------------------------------------------------

# apply NMF MS2 only
H_ms1 <- ngmcas_res$A
rownames(H_ms1) <- NULL
ngmcas_res_ms2 <- nGMCAs(
  X.m = ms2_mixedmat,
  rank = rank,
  maximumIteration = 200,
  maxFBIteration = 100,
  toleranceFB = 1e-05,
  initialization_method = "subSample",
  H_sub = H_ms1,
  errors_print = FALSE,
  method = "svsd"
)

S_ms2 <- ngmcas_res_ms2$S
pure_rt_ms2 <- melt(ngmcas_res_ms2$A) %>% as.data.table()
setnames(pure_rt_ms2, c("rank", "scan_norm", "value"))
pure_mz_ms2 <- melt(S_ms2) %>% as.data.table()
setnames(pure_mz_ms2, c("xic_label", "rank", "value"))
pure_mz_ms2 <- merge(ms2_infos, pure_mz_ms2, by = "xic_label", allow.cartesian=TRUE)

# matching
{
  NewDB <- readRDS("~/1workflow/DIA_NMF_workflow/Data/NewDB.rds")
  
  subset_rt <- 297 + c(-10, +10)
  exist_comp <- input_targ[rt_sec %between% (subset_rt)]
  ms2_mz_tolerance <- 0.007
  
  res <- lapply(5:length(unique(pure_mz_ms2$rank)), function(s){
    # s= 5
    print(paste("source", s))
    ms2_spectrum <- pure_mz_ms2[rank == s & msLevel == 2]
    
    if( nrow(ms2_spectrum) > 1 ){
      colnames(ms2_spectrum)[c(8,10)] <- c("mz_value", "intensity")
      ms2_spectrum <- ms2_spectrum[ms2_spectrum$intensity > 0, ]
      ms2_spectrum$intensity <- ms2_spectrum$intensity / max(ms2_spectrum$intensity)
      ms2_spectrum$mz_value <- gsub("\\.\\d$", "",  ms2_spectrum$mz_value) # Remove trailing '.X' where X is any digit
      ms2_spectrum$mz_value <- as.numeric( ms2_spectrum$mz_value)
      ms2_spectrum <- as.data.frame(ms2_spectrum)
      
      res2 <- lapply(1:10, function(iso){
        # iso = 1
        print(iso)
        target_iso <- unique(ms2_spectrum$isolationWindowTargetMz)[iso]
        if( !is.na(target_iso) ){
          mz_limits <- ms2_spectrum[ ms2_spectrum$isolationWindowTargetMz == target_iso,  ][1, c("lowMZ", "highMZ") ]
          ms2_spectrum_specific <- ms2_spectrum[ ms2_spectrum$isolationWindowTargetMz == target_iso,  ]
          
          if( nrow(ms2_spectrum_specific) > 1){
            c <- ress2[ ress2$source == s, ]$name
            mz <- exist_comp[Compound == c, ]$mz_pos
            ms2_scores <- match_pure_scores2(polarity = "POS", mz_prec = mz, data_base = NewDB,
                                             measured_spectra = ms2_spectrum_specific, mz_tol = ms2_mz_tolerance, mz_limits = mz_limits)
            if( nrow(ms2_scores) == 0 ){
              ms2_score <- NULL
            }else{
              lib_indx <- which.max(ms2_scores$total_score)
              ms2_score <- ms2_scores[lib_indx, ]
            }
            
            d.plot <- file.path('~/DIA_NMF_R_package/Results2/ms1/')
            if (!dir.exists(d.plot)) { dir.create(d.plot, recursive = TRUE) }
            res_sub1 <- ms2_score %>% 
              arrange(desc(total_score)) %>%
              slice_head(n = 1)
            
            lib_spectrum <- NewDB$spectra[[res_sub1$ref_spectrum_index]]@spectrum
            colnames(lib_spectrum) <- c("mz_value", "intensity") 
            lib_spectrum <- as.data.frame(lib_spectrum)
            lib_spectrum$intensity <- lib_spectrum$intensity / max(lib_spectrum$intensity)
            lib_spectrum <- lib_spectrum[ lib_spectrum$mz_value %between% mz_limits, ]
            
            p <- ggplot2::ggplot() +
              geom_linerange(data = ms2_spectrum_specific, 
                             aes(x = mz_value, ymin = 0, ymax = intensity), color = 'blue') +
              geom_linerange(data = lib_spectrum, aes(x = mz_value, ymin = -intensity, ymax = 0), color = 'red') +
              labs(x = "mz",
                   y = "intensity",
                   caption = paste('DP =', res_sub1$dot_prod, 
                                   '| IDP =', res_sub1$rev_prod, 
                                   '| FP =', res_sub1$presence, 
                                   '| Name =', res_sub1$name.))
            ggsave(paste0(d.plot, '/', 'source_', s, '_isoWind_', iso, '_ms2.png'), p, w=15, h = 10, dpi = 300)
            
            return(res_sub1)
          } else {
            return(NULL)
          }
        }
      })
    }
  })
  return(res)
  ress <- res
  names(ress) <- paste0("source_", seq(1,rank))
  ress <- do.call(rbind, ress)
}


#---------------------------------------------------------------------------------

image(ngmcas_res$A, col = scales::viridis_pal(option = "A")(30))
image(ngmcas_res$S, col = scales::viridis_pal(option = "A")(30))

# }

### Extract raw spectra ----------------------------
### Assemble EIC ----------------------------
### NMF ----------------------------
### Export results ----------------------------
res_bench <- lapply(c(10, 100, 200, 300, 500, 600, 700, 800, 1000, 1500, 2000, 5000), function(x) {
  in_mat <- mixedmat[sample(seq_len(nrow(mixedmat)), size = x, replace = T),]
  print(x)
  meantime <- microbenchmark::microbenchmark(
    "svds" = {
      nGMCAs(
        X.m = in_mat,
        rank = rank,
        maximumIteration = 5,
        maxFBIteration = 20,
        toleranceFB = 1e-05,
        initialization_method = "nndsvd",
        errors_print = FALSE,
        method = "svsd"
      )
    },
    times = 3
  )$time[3]
  data.table("n" = x, "time" = meantime)
}) %>%
  rbindlist()
plot(res_bench$n, res_bench$time / 1e9)

-