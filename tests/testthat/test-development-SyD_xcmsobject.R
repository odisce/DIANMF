testthat::skip("development script to implement xcms sequence objects")

## Improve XCMS parameters
## Run on Barbier data

library(MSnbase)
library(data.table)
input_dir <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"
input_target <- "/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Compound_data.csv"


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
for (i in seq_len(ms1_peaks[, .N])) {
  #  i <- 1
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
  
  ## build xics from peak list
  xic_dt <- peaks_i[, {
    rtrange <- c(rtmin, rtmax)
    mzrange <- c(mzmin, mzmax)
    raw_dt[mz %between% mzrange & rtime %between% rtrange, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, collisionEnergy)]
  }, by = .(peakid, peakfull)]
  xic_dt[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
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
    plot_spectrum <- ggplot(peaks_i, aes(mz, into, group = peakid)) +
      geom_linerange(aes(color = peakfull, ymin = 0, ymax = into)) +
      geom_linerange(data = peak_i, color = "black", aes(color = peakfull, ymin = 0, ymax = into)) +
      facet_grid("1" ~ .) +
      theme_bw() +
      guides(color = "none") +
      scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    plot_peak_map
    plot_xics
    plot_spectrum
  }

  ## build matrix (only on xcms peaks for now)
  ms1_mixeddt <- dcast(xic_dt[msLevel == 1, ], peakid ~ rtime, value.var = "intensity", fun.aggregate = max, fill = 0)
  ms1_mixedmat <- as.matrix(ms1_mixeddt[, -1])
  ## NMF
  rank <- 20
  ngmcas_res <- nGMCAs(
    X.m = ms1_mixedmat,
    rank = rank,
    maximumIteration = 5,
    maxFBIteration = 20,
    toleranceFB = 1e-05,
    initialization_method = "nndsvd",
    errors_print = FALSE
  )
  {
    pure_rt <- melt(ngmcas_res$A) %>% as.data.table()
    pure_mz <- melt(ngmcas_res$S) %>% as.data.table()
    ms1_mixeddt[, mz_id := seq_len(.N)]
    mz_dic <- merge(
      ms1_mixeddt[, .(peakid, mz_id)], peaks_i[, .(peakid, mz)], by = "peakid"
    )
    pure_mz <- merge(
      pure_mz[, .(mz_id = Var1, rank = Var2, value)],
      mz_dic,
      by = "mz_id"
    )

    ggpubr::ggarrange(
      ggpubr::ggarrange(
        plot_xics,
        plot_spectrum,
        align = "hv",
        ncol = 2
      ),
      ggpubr::ggarrange(
        ggplot(pure_rt, aes(Var2, value, group = as.factor(Var1))) +
          geom_line() +
          facet_grid(Var1 ~ ., scales = "free_y") +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
        ggplot(pure_mz, aes(mz, value, group = as.factor(rank))) +
          geom_linerange(aes(ymin = 0, ymax = value)) +
          facet_grid(rank ~ .) +
          theme_bw() +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
        align = "hv",
        ncol = 2
      ),
      nrow = 2,
      heights = c(1, 2)
    )
  }
  
  image(ngmcas_res$A, col = scales::viridis_pal(option = "A")(30))
  image(ngmcas_res$S, col = scales::viridis_pal(option = "A")(30))
  
}

### Extract raw spectra ----------------------------
### Assemble EIC ----------------------------
### NMF ----------------------------
### Export results ----------------------------
