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
input_dir <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
input_target <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Compound_data.csv"
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
register(MulticoreParam(10))
bpparam()

## Detect LCpeaks ----------------------------
temp_saveL <- T
# subset_rt <- 348 + c(-20, +20) # mz = 611.16
# subet_rt <- 498 + c(-20, +20) # 498 mz=500.303
subset_rt <- 49.8 + c(-20, +20) # mz = 138.055
input_targ[rt_sec %between% (subset_rt)]

target_mz <- 138.055
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
  useapexL <- T
  peak_i <- ms1_peaks[i, ]
  if (useapexL) {
    ## use apex window
    scan_rt <- 10
    rt_range <- peak_i[, (rt + c(-scan_rt, +scan_rt))]
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
      geom_pointrange(data = peak_i, aes(x = rt, xmin = rtmin, xmax = rtmax), color = "red") +
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
  raw_dt <- xcms::filterRt(xcms_obj, rt_range) %>%
    xcms::filterFile(., 1) %>%
    # xcms::filterMsLevel(., 1) %>%
    get_spectra_values()
  ## Normalize time
  time_dic <- raw_dt[, .(rtime, msLevel, isolationWindowTargetMz)] %>% unique()
  time_dic[order(rtime), scan_norm := seq_len(.N), by = .(msLevel, isolationWindowTargetMz)]
  raw_dt <- merge(
    raw_dt,
    time_dic,
    by = c("rtime", "msLevel", "isolationWindowTargetMz")
  )
  
  ## build MS1 xics from peak list
  xic_dt <- peaks_i[, {
    mzrange <- c(mzmin, mzmax)
    raw_dt[msLevel == 1 & mz %between% mzrange, .(mz, scan_norm, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy)]
  }, by = .(peakid, peakfull)]
  xic_dt[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
 
  ## build MS2 xics from raw data
  ### Combine all MS2 spectra
  MS2_peaklist <- xcms::filterRt(xcms_obj, rt_range) %>%
    xcms::filterFile(., 1) %>%
    xcms::filterMsLevel(., 2) %>%
    xcms::spectra() %>%
    Spectra::combineSpectra(ppm = 10, tolerance = 0.01, minProp = 0.1) %>%
    Spectra::asDataFrame() %>%
    data.table::as.data.table() %>%
    {.[, .(msLevel, mz, rtime, intensity)]}
  MS2_peaklist[, xic_label := paste0(seq_len(.N), "-", 2)]
  ### Extract XICs
  xic_dt_ms2 <- MS2_peaklist[, {
    mzrange <- PpmRange(mz, 30)
    raw_dt[msLevel == 2 & mz %between% mzrange, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy)]
  }, by = .(xic_label)]
  xic_dt_ms2 <- xic_dt_ms2[xic_label %in% xic_dt_ms2[, .N, by = xic_label][N >= 4, xic_label]]
  xic_dt_ms2 <- merge(
    xic_dt_ms2,
    time_dic,
    by = c("rtime", "msLevel", "isolationWindowTargetMz")
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
    plot_peak_map <- ggplot(xic_dt, aes(rtime, mz, group = peakid)) +
      geom_line(aes(color = peakfull)) +
      geom_line(data = xic_dt[peakid == peak_i[, peakid], ], color = "black") +
      facet_grid(msLevel ~ .) +
      theme_bw()
    plot_xics <- ggplot(xic_dt, aes(rtime, intensity, group = peakid)) +
      geom_line(aes(color = peakfull)) +
      geom_line(data = xic_dt[peakid == peak_i[, peakid], ], color = "black") +
      facet_grid(msLevel ~ .) +
      theme_bw() +
      guides(color = "none") +
      scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    plot_xics_ms2 <- ggplot(xic_dt_ms2, aes(rtime, intensity, group = xic_label)) +
      geom_line(aes(color = as.factor(isolationWindowTargetMz))) +
      facet_grid(isolationWindowTargetMz ~ .) +
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
    plot_peak_map
    plot_spectrum
  }

  ## build matrix (only on xcms peaks for now)
  MS1MS2_L <- TRUE
  ms1_mixeddt <- dcast(xic_dt[msLevel == 1, ], xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
  ms1_mixedmat <- as.matrix(ms1_mixeddt, rownames = TRUE)
  ms1_infos <- xic_dt[msLevel == 1, ][, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]
  if (MS1MS2_L) {
    ms2_mixeddt <- dcast(xic_dt_ms2, xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
    ms2_mixedmat <- as.matrix(ms2_mixeddt, rownames = TRUE)
    mixedmat <- rbind(
      ms1_mixedmat,
      ms2_mixedmat
    )
    ms_infos <- rbind(
      ms1_infos,
      ms2_infos
    )
    dim(mixedmat)
    ms2_infos <- xic_dt_ms2[, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]
  } else {
    mixedmat <- ms1_mixedmat
    ms_infos <- ms1_infos
  }
  REMOVE_COLINEARITY_L <- FALSE

  ## NMF
  rank <- 20
  ngmcas_res <- nGMCAs(
    X.m = mixedmat,
    rank = rank,
    maximumIteration = 5,
    maxFBIteration = 20,
    toleranceFB = 1e-05,
    initialization_method = "nndsvd",
    errors_print = FALSE,
    method = "svsd"
  ) 
  {
    pure_rt <- melt(ngmcas_res$A) %>% as.data.table()
    setnames(pure_rt, c("rank", "scan_norm", "value"))
    pure_mz <- melt(ngmcas_res$S) %>% as.data.table()
    setnames(pure_mz, c("xic_label", "rank", "value"))
    pure_mz <- merge(ms_infos, pure_mz, by = "xic_label")
    

    ggpubr::ggarrange(
      ggpubr::ggarrange(
        plot_xics,
        plot_spectrum,
        align = "hv",
        ncol = 2
      ),
      ggpubr::ggarrange(
        ggplot(pure_rt, aes(scan_norm, value, group = as.factor(rank))) +
          geom_line() +
          facet_grid(rank ~ ., scales = "free_y") +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
        ggplot() +
          geom_hline(yintercept = 0) +
          geom_linerange(data = pure_mz[msLevel == 1, ], aes(mz, value, ymin = 0, ymax = value)) +
          geom_linerange(data = pure_mz[msLevel == 2, ], aes(mz, -value, ymin = 0, ymax = -value), color = "red") +
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
