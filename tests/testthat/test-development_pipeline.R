library(MSnbase)
library(data.table)
library(BiocParallel)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)

devtools::load_all()
input_dir <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzml/"

# Step 1:
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

# register(SerialParam())
register(MulticoreParam(3))
bpparam()

## Detect LCpeaks ---------------------------- # On the complete file
temp_saveL <- T
if (temp_saveL) {
  save_path <- "./temp2/data/"
  dir.create(save_path, recursive = TRUE)
  xcms_obj_path <- file.path(save_path, "xcms_obj.rds")
  if (file.exists(xcms_obj_path)) {
    xcms_obj <- readRDS(xcms_obj_path)
  } else {
    xcms_obj <- detect_xcms_peaks(
      sequence_table = mzml_dt[, lapply(.SD, head, 1), by = .(class)],
      params = params_ls
      # rt_range = subset_rt
    )
    saveRDS(xcms_obj, file = xcms_obj_path)
  }
}


## extract all detected MS1 peaks, from all samples ----------------------------
ms1_peaks <- extract_xcms_peaks(xcms_obj)
ms1_peaks[, iteration := as.integer(NA) ]

# extract aligned features; we will work on them
ms1_features <- extract_xcms_features(xcms_obj, orderL = T, orderL_sample = "20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML")
sample_idx <- 1
sample_name <- "20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML"
ms1_features_sample <- ms1_features[, .SD, .SDcols = c("featureid", "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "peakidx", sample_name)]
colnames(ms1_features_sample) <- c("featureid", "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "peakidx", sample_name)
ms1_features[, iteration := as.integer(NA) ]

min_rt <- min(ms1_peaks[sample == sample_idx, ]$rtmin)
max_rt <- max(ms1_peaks[sample == sample_idx, ]$rtmax)

features.l <- list()
feature_idx <- 27
k <- 1

while( feature_idx <= nrow(ms1_features) ) {
  
  print(paste("feature index:", feature_idx, '------ k:', k))
  peaks_idxs <- unlist(ms1_features[feature_idx, ]$peakidx)
  # peaks_idxs <- paste0("CP", sprintf("%06d", peaks_idxs))
  ms1_peaks_sub <- ms1_peaks[peaks_idxs, ]
  ms1_peak <- ms1_peaks_sub[sample == sample_idx, ]
  ms1_peak <- ms1_peak[which.max(ms1_peak$into), ]
  
  if( nrow(ms1_peak) != 0 ){
    target_mz <- ms1_peak$mz
    useapexL <- T # -------------------------------------------------------------- param
    peak_i <- ms1_peak
    if(useapexL) {
      ## use apex window
      scan_rt <- 10
      rt_range <- peak_i[, (rt + c(-scan_rt-1, +scan_rt+1))]
      rt_range[1] <- max(min_rt, rt_range[1])  # to avoid rt < 0
      rt_range[2] <- min(max_rt, rt_range[2])  # to avoid rt out of range
      peaks_i <- ms1_peaks[sample == sample_idx & rtmin <= rt_range[2] & rtmax >= rt_range[1], ]
    } else {
      ## full peak range
      rt_range <- peak_i[, c(rtmin, rtmax)]
      peaks_i <- ms1_peaks[sample == sample_idx & rtmin <= peak_i[, rtmax] & rtmax >= peak_i[, rtmin], ]
    }
    
    ## flag peaks partially, fully or apex inside the window
    peaks_i[, peakfull := ifelse(
      (rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
      ifelse(rt %between% rt_range, "apex", "partial")
    )]
    { ## plot peaks in range
      ggplot(peaks_i, aes(rt, mz, group = peakid)) +
        geom_pointrange(aes(x = rt, xmin = rtmin, xmax = rtmax, color = peakfull)) +
        geom_vline(xintercept = rt_range, linetype = 2, color = "red") +
        theme_bw() +
        labs(
          title = paste0("MS1 peaks in rt range: [", round(rt_range[1],4), ',', round(rt_range[2], 4), ']'),
          caption = sprintf(
            "Total MS1 peaks: %i",
            peaks_i[, .N] )
        )
    }

    # extract data
    raw_dt <- xcms::filterRt(xcms_obj, rt_range) %>%
      xcms::filterFile(., sample_idx) %>%
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
    
    ## build MS2 xics from raw data
    MS2_ISOEACHL <- T  # from each isolation window separately # ----------------- param
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
    
    ### Extract MS2 XICs
    xic_dt_ms2 <- MS2_peaklist[, {
      mzrange <- PpmRange(mz, 7)
      raw_dt[msLevel == 2 & mz %between% mzrange, .(mz, rtime, intensity, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz,
                                                    collisionEnergy)]
    }, by = .(xic_label)]
    xic_dt_ms2 <- xic_dt_ms2[xic_label %in% xic_dt_ms2[, .N, by = xic_label][N >= 4, xic_label]]
    xic_dt_ms2 <- merge(
      xic_dt_ms2,
      time_dic,
      by = c("rtime", "msLevel", "isolationWindowTargetMz"),
      allow.cartesian = T
    )
    xic_dt_ms2[, xic_label := paste0(xic_label, "-", isolationWindowTargetMz)]
    
    ## build matrix
    MS1MS2_L <- F #---------------------------------------------------------------param
    # ms1 matrix from xcms eics
    ms1_mixeddt <- dcast(xic_dt_ms1[msLevel == 1, ], xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
    ms1_mixedmat <- ms1_mixeddt <- as.matrix(ms1_mixeddt, rownames = TRUE)
    row_filter_ms1 <- apply(ms1_mixedmat, 1, has_four_consecutive_non_zero)
    ms1_mixedmat <- ms1_mixedmat[row_filter_ms1, , drop = FALSE]
    ms1_mixedmat_deleted <- ms1_mixeddt[!row_filter_ms1, , drop = FALSE]
    
    ms1_infos <- xic_dt_ms1[msLevel == 1, ][, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz)]
    ms1_infos$isolationWindowLowerMz <- NA
    ms1_infos$isolationWindowUpperMz <- NA
    
    ms2_mixeddt <- dcast(xic_dt_ms2, xic_label ~ scan_norm, value.var = "intensity", fun.aggregate = max, fill = 0)
    ms2_mixedmat <- as.matrix(ms2_mixeddt, rownames = TRUE)
    row_filter_ms2 <- apply(ms2_mixedmat, 1, has_four_consecutive_non_zero)
    ms2_mixedmat <- ms2_mixedmat[row_filter_ms2, , drop = FALSE]
    ms2_infos <- xic_dt_ms2[, .(mz = median(mz)), by = .(xic_label, msLevel, isolationWindowTargetMz, isolationWindowLowerMz, isolationWindowUpperMz)]
    
    if(MS1MS2_L) {
      mixedmat <- rbind(
        ms1_mixedmat,
        ms2_mixedmat )
      
      ms_infos <- rbind(
        ms1_infos,
        ms2_infos )  }
    
    ## NMF
    if(MS1MS2_L == F){ # MS1 then MS2 separately
      
      ### NMF on MS1
      rank <- min(15, ncol(ms1_mixedmat))  # to be changed
      ngmcas_res_ms1 <- nGMCAs(
        X.m = ms1_mixedmat,
        rank = rank,
        maximumIteration = 200,
        maxFBIteration = 100,
        toleranceFB = 1e-05,
        initialization_method = "nndsvd",
        errors_print = FALSE,
        method = "svsd"
      )
      
      pure_rt_ms1 <- melt(ngmcas_res_ms1$A) %>% as.data.table()
      setnames(pure_rt_ms1, c("rank", "scan_norm", "value"))
      S <- ngmcas_res_ms1$S
      pure_mz_ms1 <- melt(S) %>% as.data.table()
      setnames(pure_mz_ms1, c("xic_label", "rank", "value"))
      pure_mz_ms1 <- merge(ms1_infos, pure_mz_ms1, by = "xic_label")
      
      ### NMF on MS2
      H_ms1 <- ngmcas_res_ms1$A
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
      
    } else { # MS1 and MS2 combined
      
      ### NMF
      rank <- min(20, ncol(ms1_mixedmat)) # to be changed
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
      
      #### extract the MS1 data
      S_ms1 <- ngmcas_res$S[1:dim(ms1_mixedmat)[1], ]
      pure_rt_ms1 <- melt(ngmcas_res$A) %>% as.data.table()
      setnames(pure_rt_ms1, c("rank", "scan_norm", "value"))
      pure_mz_ms1 <- melt(S_ms1) %>% as.data.table()
      setnames(pure_mz_ms1, c("xic_label", "rank", "value"))
      pure_mz_ms1 <- merge(ms1_infos, pure_mz_ms1, by = "xic_label")
      
      #### extract the MS2 data
      S_ms2 <- ngmcas_res$S[(dim(ms1_mixedmat)[1]+1):( (dim(ms1_mixedmat)[1]) + dim(ms2_mixedmat)[1]), ]
      pure_rt_ms2 <- melt(ngmcas_res$A) %>% as.data.table()
      setnames(pure_rt_ms2, c("rank", "scan_norm", "value"))
      pure_mz_ms2 <- melt(S_ms2) %>% as.data.table()
      setnames(pure_mz_ms2, c("xic_label", "rank", "value"))
      pure_mz_ms2 <- merge(ms2_infos, pure_mz_ms2, by = "xic_label", allow.cartesian=TRUE)
      
    }
    
    # Plot extracted sources
    # {
    # plot_xics_ms1 <- ggplot(xic_dt_ms1, aes(rtime, intensity, group = peakid)) +
    #   geom_line(aes(color = peakfull)) +
    #   geom_point(aes(color = peakfull)) +
    #   theme_bw() +
    #   guides(color = "none") +
    #   scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    # 
    #   plot_xics_ms2 <- ggplot(xic_dt_ms2, aes(rtime, intensity, group = xic_label)) +
    #     geom_line(aes(color = xic_label)) +
    #     geom_point(aes(color = xic_label)) +
    #     theme_bw() +
    #     guides(color = "none") +
    #     scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    # 
    #   plot_spectrum_raw <- ggplot(peaks_i, aes(mz, into, group = peakid)) +
    #     geom_linerange(aes(color = peakfull, ymin = 0, ymax = into)) +
    #     geom_linerange(data = peak_i, color = "black", aes(color = peakfull, ymin = 0, ymax = into)) +
    #     theme_bw() +
    #     guides(color = "none") +
    #     scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    #   # plot mixed MS1 (after filtering)
    #   ms1_mixed_spectra <- ms1_mixedmat
    #   rownames(ms1_mixed_spectra) <- sub("-.*", "",  rownames(ms1_mixed_spectra))
    #   rownames(ms1_mixed_spectra) <- ms1_peaks[peakid %in%  rownames(ms1_mixed_spectra), ]$mz
    #   ms1_mixed_spectra <- reshape2::melt(ms1_mixed_spectra)
    #   colnames(ms1_mixed_spectra) <- c('mz', 'rt', 'into')
    #   ms1_mixed_spectra <- ms1_mixed_spectra[ms1_mixed_spectra$into > 0, ]
    #   plot_spectrum <- ggplot() +
    #     geom_linerange(data = ms1_mixed_spectra, aes(x = mz, ymin = 0, ymax = into)) +
    #     theme_bw(base_size = 14)
    # 
    #   p1 <- ggpubr::ggarrange(
    #     ggpubr::ggarrange(
    #       plot_xics_ms1,
    #       plot_xics_ms2,
    #       plot_spectrum_raw,
    #       plot_spectrum,
    #       align = "hv",
    #       ncol = 4
    #     ),
    #     ggpubr::ggarrange(
    #       ggplot(pure_rt_ms1, aes(scan_norm, value, group = as.factor(rank))) +
    #         geom_line() +
    #         facet_grid(rank ~ ., scales = "free_y") +
    #         theme_bw() +
    #         scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
    #       ggplot() +
    #         geom_hline(yintercept = 0) +
    #         geom_linerange(data = pure_mz_ms1[msLevel == 1, ], aes(mz, value, ymin = 0, ymax = value)) +
    #         geom_linerange(data = pure_mz_ms2[msLevel == 2, ], aes(mz, -value, ymin = 0, ymax = -value), color = "red") +
    #         facet_grid(rank ~ .) +
    #         theme_bw() +
    #         scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    #         guides(color = "none"),
    #       align = "hv",
    #       ncol = 2
    #     ),
    #     nrow = 2,
    #     heights = c(1, 2)
    #   )
    #   # p1
    #   # ggplot2::ggsave(paste0('~/DIA_NMF_R_package/Results3/', k, '_sources.png'), p1, w=10, h = 10, dpi = 300)
    # }

    ## quantifying every xcms peak by its pure spectra and delete noisy sources
    peaks_xcms <- peaks_i  
    peaks_xcms$peakid <- paste0(peaks_xcms$peakid, "-1")
    colnames(peaks_xcms)[1] <- "xic_label"
    
    setDT(peaks_xcms)
    setDT(pure_mz_ms1)
    merge_data <- peaks_xcms[pure_mz_ms1, on = "xic_label"]
    merge_data <- merge_data[, c("xic_label", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "peakfull" , "i.mz", "rank", "value")]
    
    merge_data <- merge_data %>%
      group_by(xic_label) %>%
      mutate(contribution = value / sum(value)) %>%
      mutate(contribution = replace(contribution, is.nan(contribution), 0)) %>%
      ungroup()
    
    setDT(merge_data)
    contribution_matrix <- dcast(
      merge_data, 
      rank ~ xic_label, 
      value.var = "contribution", 
      fill = 0
    )
    contribution_matrix <- as.matrix(contribution_matrix, rownames = TRUE)
    contribution_matrix <- t(contribution_matrix)
    
    # find sources who have max contribution higher than 0.6
    good_sources <- which( unname(apply(contribution_matrix, 2, max)) >=0.6 )
    
    ### find where every peak contribute the most
    peaks_sources_df <- data.frame(
      xic_label = rownames(contribution_matrix),
      source = colnames(contribution_matrix)[apply(contribution_matrix, 1, which.max)]
    )
    
    info <- merge(peaks_sources_df, peaks_xcms, by = "xic_label")
    peaks_p <- ggplot(info, aes(rt, mz, group = xic_label)) +
      geom_pointrange(aes(x = rt, xmin = rtmin, xmax = rtmax, color = source)) +
      geom_vline(xintercept = rt_range, linetype = 2, color = "red") +
      # facet_grid(source ~ .) +
      theme_bw()
    # peaks_p <- plotly::ggplotly(peaks_p)
    # ggsave(paste0('~/DIA_NMF_R_package/Results3/', k, '_peaks.png'), peaks_p, w=10, h = 10, dpi = 300)
    
    # get the real rt 
    peak_counts <- xic_dt_ms1[, .N, by = peakid]
    peak_id <- which(peak_counts$N == ncol(ms1_mixedmat))[1]
    peak_id <- peak_counts[peak_id, ]$peakid
    real_rt <-  xic_dt_ms1[ peakid == peak_id, ]$rtime 
    chroms <- lapply(1:rank, function(eic){
      ints <- pure_rt_ms1[rank == eic, ]$value
      rt <- real_rt
      ch <- Chromatogram(rtime = rt, ints)
    })
    
    chrs <- MChromatograms(chroms, nrow = rank)
    detected_peaks <- findChromPeaks(object = chrs,
                                     param =  CentWaveParam(
                                       ppm = 0,
                                       peakwidth = c(4, 10),
                                       snthresh = 0,
                                       mzCenterFun = "wMeanApex3",
                                       integrate = 2,
                                       mzdiff = -0.0003,
                                       noise = 0,
                                       firstBaselineCheck = FALSE),
                                     msLevel = 1L)
    peaks <- chromPeaks(detected_peaks)
    peaks <- as.data.frame(peaks)
    peaks <- peaks %>%  # filter sources of more than 2 peaks, and of negative sn
      group_by(row) %>%  
      filter(n() <= 2 & all(sn >= 0)) %>%
      ungroup()
    peaks <- setDT(peaks)
    peaks <- peaks[ abs(rt - rt_range[1]) >= 3 & abs(rt - rt_range[2]) >= 3, ] # filter peaks of rt very close to one edge of rt_range
    
    # filter peaks, which doesn't contains peaks of high contribution (lower 0.6)
    peaks <- peaks[ row %in% good_sources, ]
    
    # filter peaks 
    if( nrow(peaks) > 0 ){
      res <- lapply(1:nrow(peaks), function(r){
        s <- peaks[r,  ]$row
        info_sub <- info[info$source == s, ]
        info_sub <- info_sub[ info_sub$rt %between% c(peaks[r, ]$rtmin, peaks[r, ]$rtmax),  ]
      })
      info_new <- do.call(rbind, res)
      
      if( nrow(info_new) > 0 ){
        peaks_p <- ggplot(info_new, aes(rt, mz, group = xic_label)) +
          geom_pointrange(aes(x = rt, xmin = rtmin, xmax = rtmax, color = source)) +
          geom_vline(xintercept = rt_range, linetype = 2, color = "red") +
          # facet_grid(source ~ .) + 
          theme_bw()
        # peaks_p <- plotly::ggplotly(peaks_p)
        # ggsave(paste0('~/DIA_NMF_R_package/Results3/', k, '_selected_peaks.png'), peaks_p, w=10, h = 10, dpi = 300)
        
        ## Extract MS1 and MS2 pure sources
        good_sources <- as.numeric(unique(info_new$source))
        pure_mz_ms1 <- pure_mz_ms1[ pure_mz_ms1$rank %in% good_sources, ]
        pure_mz_ms1 <- pure_mz_ms1[ value > 0, ]  # filter zero intensity fragments
        pure_mz_ms2 <- pure_mz_ms2[ pure_mz_ms2$rank %in% good_sources, ]
        pure_mz_ms2 <- pure_mz_ms2[ value > 0, ]
        pure_rt_ms1 <- pure_rt_ms1[ pure_rt_ms1$rank %in% good_sources, ]
        pure_rt_ms2 <- pure_rt_ms2[ pure_rt_ms2$rank %in% good_sources, ]
        
        # {
        #   p2 <- ggpubr::ggarrange(
        #     ggpubr::ggarrange(
        #       plot_xics_ms1,
        #       plot_xics_ms2,
        #       plot_spectrum,
        #       align = "hv",
        #       ncol = 3
        #     ),
        #     ggpubr::ggarrange(
        #       ggplot(pure_rt_ms1, aes(scan_norm, value, group = as.factor(rank))) +
        #         geom_line() +
        #         facet_grid(rank ~ ., scales = "free_y") +
        #         theme_bw() +
        #         scale_y_continuous(labels = function(x) format(x, scientific = TRUE)),
        #       ggplot() +
        #         geom_hline(yintercept = 0) +
        #         geom_linerange(data = pure_mz_ms1[msLevel == 1, ], aes(mz, value, ymin = 0, ymax = value)) +
        #         geom_linerange(data = pure_mz_ms2[msLevel == 2, ], aes(mz, -value, ymin = 0, ymax = -value), color = "red") +
        #         facet_grid(rank ~ .) +
        #         theme_bw() +
        #         scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
        #         guides(color = "none"),
        #       align = "hv",
        #       ncol = 2
        #     ),
        #     nrow = 2,
        #     heights = c(1, 2)
        #   )
        #   # p2
        #   # ggsave(paste0('~/DIA_NMF_R_package/Results3/', k, '_selected_sources.png'), p2, w=10, h = 10, dpi = 300)
        # }

        feature_sub.l <- list(
          'feature' = ms1_features[feature_idx, ],
          'peak' = ms1_peak,
          'MS1_mixed_mat' = ms1_mixedmat,
          'MS2_mixed_mat' = ms2_mixedmat,
          'ngmca_res' = if(MS1MS2_L) ngmcas_res else list(ngmcas_res_ms1, ngmcas_res_ms2),
          'ms1_pure_spectra' = pure_mz_ms1,  
          'ms2_pure_spectra' = pure_mz_ms2,
          'MS1_pure_elution_profiles' = pure_rt_ms1,   
          'MS2_pure_elution_profiles' = pure_rt_ms1,
          'xcms_assigned_sources' = info_new,
          'ms1_info' = list('ms1_infos' = ms1_infos, 'xic_dt_ms1' = xic_dt_ms1),
          'ms2_info' = list('ms2_infos' = ms2_infos, 'xic_dt_ms2' = xic_dt_ms2)  
          )
        
        features.l[[k]] <- feature_sub.l
        
        idx_deleted <- paste0(ms1_peaks$peakid, "-1") %in% rownames(ms1_mixedmat_deleted)  # eics, which doesn't have at least 4 consecutive non-zero values
        
        ## all peaks fully included in the rt range and mostly included with their apex in the rt range should be not processed again
        ## since, these peaks are well detected peaks (so well extracted in the pure sources) or, noises (not included in the pure sources, so we should not return for this rt range for them)
        min_distance <- 5 #--------------------------------------------------------- param
        filtered_peaks <- peaks_xcms %>%
          filter(
            peakfull == "full" | 
              (peakfull == "apex" & (abs(rt-rt_range[1]) >= min_distance) & abs(rt-rt_range[2]) >= min_distance) )
        idx_sub <- paste0(ms1_peaks$peakid, "-1") %in% filtered_peaks$xic_label
        
        combined_idx <- c(which(idx_deleted), which(idx_sub))
        idx <- unique(combined_idx)
        ms1_peaks$iteration[idx] <- ifelse(
          is.na(ms1_peaks$iteration[idx]), 
          k, 
          paste0(ms1_peaks$iteration[idx], ",", k)  )
      }
    }
    
    if( nrow(peaks) == 0 | (exists("info_new") & nrow(info_new) == 0) ) {
      ## update the peaks df
      idx <- paste0(ms1_peaks$peakid, "-1") %in% peaks_xcms$xic_label  # all detected peaks
      idx <- which(idx)
      ms1_peaks$iteration[idx] <- ifelse(
        is.na(ms1_peaks$iteration[idx]),
        0,
        paste0(ms1_peaks$iteration[idx], ",", 0)  ) 
    }
    
    # iterate over the next non-processed feature, but update it first
    peaks_removed <- idx
    ms1_features[, iteration := as.character(iteration)]
    ms1_features[sapply(peakidx, function(p) any(unlist(p) %in% peaks_removed)), 
                 iteration := ifelse(is.na(iteration), as.character(k), paste0(iteration, ",", k))]
    
    if( nrow(peaks) > 0 & exists("info_new") & nrow(info_new) > 0 ){
      k <- k + 1
    }
    
  } else {
    ms1_features[feature_idx, iteration := 0] 
  }
  
  feature_idx <- which(is.na(ms1_features$iteration))[1]
  if( is.na(feature_idx) ){
    break  }
  
}
  
# saveRDS(ms1_peaks, "~/DIA_NMF_R_package/Results3/ms1_peaks.rds")
# saveRDS(ms1_features, "~/DIA_NMF_R_package/Results3/ms1_features.rds")
# saveRDS(features.l, "~/DIA_NMF_R_package/Results3/features_list.rds")

# saveRDS(pure_mz_ms1, "~/DIA_NMF_R_package/pure_mz_ms1.rds")


# library(dplyr)
# ms1_peaks <- ms1_peaks %>%
#   mutate(iteration = sub(",.*", "", iteration))  # Keep only the first value before the comma
# ggplot(ms1_peaks, aes(rt, mz, group = iteration)) + 
#   geom_pointrange(aes(x = rt, xmin = rtmin, xmax = rtmax, color = iteration)) +
#   theme_bw()+
#   guides(color = "none")
