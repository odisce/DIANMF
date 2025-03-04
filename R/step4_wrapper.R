# wrapper function
## xcms parameters
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



#' DIANMF wrapper function
#'
#' @param input_dir `character` mzML files directory path.
#' @param d.out  `character` output directory path.
#' @param sample_idx `numeric` index of the sample.
#' @param temp_saveL `logical`
#' @param MS2_ISOEACHL `logical`
#' @param MS1MS2_L `logical`
#' @param scan_rt_ext `numeric`
#' @param min_distance `numeric`
#' @param plots `logical`
#' @inheritParams detect_LCfeatures
#' 
#' @return `list`
#' @export
#'
#' @import xcms
#' @import MSnbase
#' @import data.table
#' @import magrittr
#' @import dplyr
DIANMF_f <- function(input_dir, d.out, 
                     sample_idx = 1, 
                     temp_saveL = T, MS2_ISOEACHL = T, MS1MS2_L = F,
                     scan_rt_ext = 10, min_distance = 5,
                     params_ls, plots = F ){
  
  mzml_dt <- prepare_mzMLfiles(input_dir)
  
  xcms_obj <- detect_LCfeatures(temp_saveL, params_ls)
  ms1_peaks <- extract_xcms_peaks(xcms_obj)
  ms1_features <- extract_xcms_features(xcms_obj)
  ms1_features[, iteration := as.integer(NA) ]
  
  min_rt <- min(ms1_peaks[sample == sample_idx, ]$rtmin)
  max_rt <- max(ms1_peaks[sample == sample_idx, ]$rtmax)
  features.l <- list()
  feature_idx <- 1
  k <- 1
  
  # start iterate over rt_ranges 
  while( feature_idx <= nrow(ms1_features) ) {
    
    print(paste("feature index:", feature_idx, '------ k:', k))
    peaks_idxs <- unlist(ms1_features[feature_idx, ]$peakidx)
    ms1_peaks_sub <- ms1_peaks[peaks_idxs, ]
    ms1_peak <- ms1_peaks_sub[sample == sample_idx, ]
    ms1_peak <- ms1_peak[which.max(ms1_peak$into), ]
    
    if( nrow(ms1_peak) != 0 ){
      target_mz <- ms1_peak$mz
      peak_i <- ms1_peak
      scan_rt <- 10
      rt_range <- peak_i[, (rt + c(-scan_rt-1, +scan_rt+1))]
      rt_range[1] <- max(min_rt, rt_range[1])  # to avoid rt < 0
      rt_range[2] <- min(max_rt, rt_range[2])  # to avoid rt out of range
      peaks_i <- ms1_peaks[sample == sample_idx & rtmin <= rt_range[2] & rtmax >= rt_range[1], ]
      
      ## flag peaks partially, fully or apex inside the window
      peaks_i[, peakfull := ifelse(
        (rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
        ifelse(rt %between% rt_range, "apex", "partial")
      )]
      
      res_general <- get_rawD_ntime(msexp = xcms_obj, rt_range)
      raw_dt <- res_general$raw_dt
      time_dic <- res_general$time_dic
      xic_dt_ms1 <- build_ms1XICS(peaks_i)
      xic_dt_ms2 <- build_ms2XICs(MS2_ISOEACHL = T, msexp = xcms_obj, time_dic)
      
      ## Generate data and matrices
      ### ms1
      ms1Data <- ms1Info(xic_dt_ms1)
      ms1_mixedmat <- ms1Data$ms1_mixedmat
      ms1_mixedmat_deleted <- ms1Data$ms1_mixedmat_deleted
      ms1_infos <- ms1Data$ms1_infos
      
      ### ms2
      ms2Data <- ms1Info(xic_dt_ms2)
      ms2_mixedmat <- ms2Data$ms2_mixedmat
      ms2_infos <- ms2Data$ms2_infos
      
      ## NMF 
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
        rank <- min(15, ncol(ms1_mixedmat))  
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
        
        pure_rt_ms1 <- reshape2::melt(ngmcas_res_ms1$A) %>% as.data.table()
        setnames(pure_rt_ms1, c("rank", "scan_norm", "value"))
        S <- ngmcas_res_ms1$S
        pure_mz_ms1 <- reshape2::melt(S) %>% as.data.table()
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
        pure_rt_ms2 <- reshape2::melt(ngmcas_res_ms2$A) %>% as.data.table()
        setnames(pure_rt_ms2, c("rank", "scan_norm", "value"))
        pure_mz_ms2 <- reshape2::melt(S_ms2) %>% as.data.table()
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
        pure_rt_ms1 <- reshape2::melt(ngmcas_res$A) %>% as.data.table()
        setnames(pure_rt_ms1, c("rank", "scan_norm", "value"))
        pure_mz_ms1 <- reshape2::melt(S_ms1) %>% as.data.table()
        setnames(pure_mz_ms1, c("xic_label", "rank", "value"))
        pure_mz_ms1 <- merge(ms1_infos, pure_mz_ms1, by = "xic_label")
        
        #### extract the MS2 data
        S_ms2 <- ngmcas_res$S[(dim(ms1_mixedmat)[1]+1):( (dim(ms1_mixedmat)[1]) + dim(ms2_mixedmat)[1]), ]
        pure_rt_ms2 <- reshape2::melt(ngmcas_res$A) %>% as.data.table()
        setnames(pure_rt_ms2, c("rank", "scan_norm", "value"))
        pure_mz_ms2 <- reshape2::melt(S_ms2) %>% as.data.table()
        setnames(pure_mz_ms2, c("xic_label", "rank", "value"))
        pure_mz_ms2 <- merge(ms2_infos, pure_mz_ms2, by = "xic_label", allow.cartesian=TRUE)
        
      }
      
      ## Quantifying every xcms peak by its pure spectra and delete noisy sources
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
        dplyr::ungroup()
      
      setDT(merge_data)
      contribution_matrix <- dcast(
        merge_data, 
        rank ~ xic_label, 
        value.var = "contribution", 
        fill = 0
      )
      contribution_matrix <- as.matrix(contribution_matrix, rownames = TRUE)
      contribution_matrix <- t(contribution_matrix)
      
      ### find sources who have max contribution higher than 0.6
      good_sources <- which( unname(apply(contribution_matrix, 2, max)) >=0.6 )
      
      ### find where every peak contribute the most
      peaks_sources_df <- data.frame(
        xic_label = rownames(contribution_matrix),
        source = colnames(contribution_matrix)[apply(contribution_matrix, 1, which.max)]
      )
      
      info <- merge(peaks_sources_df, peaks_xcms, by = "xic_label")
      
      ## Extract the good sources
      ### get the real rt 
      peak_counts <- xic_dt_ms1[, .N, by = peakid]
      peak_id <- which(peak_counts$N == ncol(ms1_mixedmat))[1]
      peak_id <- peak_counts[peak_id, ]$peakid
      real_rt <-  xic_dt_ms1[ peakid == peak_id, ]$rtime 
      chroms <- lapply(1:rank, function(eic){
        ints <- pure_rt_ms1[rank == eic, ]$value
        rt <- real_rt
        ch <- MSnbase::Chromatogram(rtime = rt, ints)
      })
      
      chrs <- MSnbase::MChromatograms(chroms, nrow = rank)
      detected_peaks <- xcms::findChromPeaks(object = chrs,
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
      peaks <- xcms::chromPeaks(detected_peaks)
      peaks <- as.data.frame(peaks)
      peaks <- peaks %>%  # filter sources of more than 2 peaks, and of negative sn
        group_by(row) %>%  
        filter(n() <= 2 & all(sn >= 0)) %>%
        ungroup()
      peaks <- setDT(peaks)
      peaks <- peaks[ abs(rt - rt_range[1]) >= 3 & abs(rt - rt_range[2]) >= 3, ] # filter peaks of rt very close to one edge of rt_range
      
      ### filter peaks, which doesn't contains peaks of high contribution (lower 0.6)
      peaks <- peaks[ row %in% good_sources, ]
      
      ### filter peaks 
      if( nrow(peaks) > 0 ){
        res <- lapply(1:nrow(peaks), function(r){
          s <- peaks[r,  ]$row
          info_sub <- info[info$source == s, ]
          info_sub <- info_sub[ info_sub$rt %between% c(peaks[r, ]$rtmin, peaks[r, ]$rtmax),  ]
        })
        info_new <- do.call(rbind, res)
        
        if( nrow(info_new) > 0 ){
          
          #### Extract MS1 and MS2 pure sources
          good_sources <- as.numeric(unique(info_new$source))
          pure_mz_ms1 <- pure_mz_ms1[ pure_mz_ms1$rank %in% good_sources, ]
          pure_mz_ms1 <- pure_mz_ms1[ value > 0, ]  # filter zero intensity fragments
          pure_mz_ms2 <- pure_mz_ms2[ pure_mz_ms2$rank %in% good_sources, ]
          pure_mz_ms2 <- pure_mz_ms2[ value > 0, ]
          pure_rt_ms1 <- pure_rt_ms1[ pure_rt_ms1$rank %in% good_sources, ]
          pure_rt_ms2 <- pure_rt_ms2[ pure_rt_ms2$rank %in% good_sources, ]
          
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
  
  # save results
  if (!is.null(d.out)) {
    if (!dir.exists(d.out)) {dir.create(d.out)} }
  saveRDS(mzml_dt, paste0(d.out, '/mzml_dt.rds'))
  saveRDS(ms1_peaks, paste0(d.out, '/ms1_peaks.rds'))
  saveRDS(ms1_features, paste0(d.out, '/ms1_features.rds'))
  saveRDS(features.l, paste0(d.out, '/features_list.rds'))
                          
  return(features.l)
}
