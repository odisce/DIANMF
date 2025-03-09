#' Process xcms object
#'
#' @inheritParams extract_xcms_peaks
#' @param d.out `logical` TRUE to save results, else FALSE
#' @param sample_idx NULL to process all files else, `numeric` index of the sample.
#' @param MS2_ISOEACHL `logical`
#' @param MS1MS2_L `logical`
#' @inheritParams nGMCAs
#' @param scan_rt_ext `numeric`
#' @param min_distance `numeric`
#'
#' @return `list`
#' @export
#' 
#' @importFrom xcms findChromPeaks chromPeaks
#' @importFrom MSnbase MChromatograms Chromatogram
#' @importFrom reshape2 melt
#' @importFrom data.table setnames setDT
#' @import magrittr
#' @import dplyr
#' @importFrom MsExperiment sampleData
#' @importFrom tools file_path_sans_ext
DIANM.f <- function(msexp,
                    d.out = FALSE,
                    sample_idx = 1, 
                    MS2_ISOEACHL = T, MS1MS2_L = F,
                    rank = 10,
                    maximumIteration = 200,
                    maxFBIteration = 100,
                    toleranceFB = 1e-05,
                    initialization_method = "nndsvd",
                    errors_print = FALSE,
                    method = "svsd",
                    scan_rt_ext = 10, min_distance = 5 ){
  
  ms1_peaks <- extract_xcms_peaks(msexp)
  ms1_features <- extract_xcms_features(msexp)
  
  file_info <- MsExperiment::sampleData(msexp)
  
  if( is.null(sample_idx) ){
    sample_idx <- file_info$InjectionOrder
  } else{
    sample_idx <- file_info$InjectionOrder[sample_idx]
  }

  res <- lapply(sample_idx, function(s_idx){
    
    print(paste("start processing sample:", s_idx))
    ms1_peaks <- ms1_peaks[, iteration := as.integer(NA) ]
    ms1_features[, iteration := as.integer(NA) ]
    
    min_rt <- min(ms1_peaks[sample == s_idx, ]$rtmin)
    max_rt <- max(ms1_peaks[sample == s_idx, ]$rtmax)
    features.l <- list()
    feature_idx <- 1
    k <- 1
    
    # start iterate over rt_ranges 
    while( feature_idx <= nrow(ms1_features) ) {
      
      print(paste("feature index:", feature_idx, '------ k:', k))
      peaks_idxs <- unlist(ms1_features[feature_idx, ]$peakidx)
      ms1_peaks_sub <- ms1_peaks[peaks_idxs, ]
      ms1_peak <- ms1_peaks_sub[sample == s_idx, ]
      ms1_peak <- ms1_peak[which.max(ms1_peak$into), ]
      
      if( nrow(ms1_peak) != 0 ){

        target_mz <- ms1_peak$mz
        peak_i <- ms1_peak
        rt_range <- peak_i[, (rt + c(-scan_rt_ext-1, +scan_rt_ext+1))]
        rt_range[1] <- max(min_rt, rt_range[1])  # to avoid rt < 0
        rt_range[2] <- min(max_rt, rt_range[2])  # to avoid rt out of range
        peaks_i <- ms1_peaks[sample == s_idx & rtmin <= rt_range[2] & rtmax >= rt_range[1], ]
        
        ## flag peaks partially, fully or apex inside the window
        peaks_i[, peakfull := ifelse(
          (rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
          ifelse(rt %between% rt_range, "apex", "partial")
        )]
        
        res_general <- get_rawD_ntime(msexp, rt_range, s_idx)
        raw_dt <- res_general$raw_dt
        time_dic <- res_general$time_dic
        xic_dt_ms1 <- build_ms1XICS(peaks_i, raw_dt)
        xic_dt_ms2 <- build_ms2XICs(msexp, raw_dt, time_dic, rt_range, MS2_ISOEACHL = T)
        
        ## Generate data and matrices
        ### ms1
        ms1Data <- ms1Info(xic_dt_ms1)
        ms1_mixedmat <- ms1Data$ms1_mixedmat
        ms1_mixedmat_deleted <- ms1Data$ms1_mixedmat_deleted
        ms1_infos <- ms1Data$ms1_infos
        
        ### ms2
        ms2Data <- ms2Info(xic_dt_ms2)
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
        contribution_matrix <- data.table::dcast(
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
        real_rt <-  time_dic[msLevel == 1, ]$rtime
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
              'MS1_mixed_mat' = ms1_mixedmat,
              'MS2_mixed_mat' = ms2_mixedmat,
              'ms1_pure_spectra' = pure_mz_ms1,  
              'ms2_pure_spectra' = pure_mz_ms2,
              'MS1_pure_elution_profiles' = pure_rt_ms1,   
              'MS2_pure_elution_profiles' = pure_rt_ms1,
              'xcms_assigned_sources' = info_new
            )
            
            features.l[[k]] <- feature_sub.l
            
            ## update the peaks df
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
            updated_data <- update_df(ms1_peaks, ms1_features, idx, iteration_number = k)
            ms1_peaks <- updated_data$ms1_peaks
            ms1_features <- updated_data$ms1_features
            
            k <- k + 1
            }else{
              
              idx <- paste0(ms1_peaks$peakid, "-1") %in% peaks_xcms$xic_label  # all detected peaks
              idx <- which(idx)
              updated_data <- update_df(ms1_peaks, ms1_features, idx, iteration_number = 0)
              ms1_peaks <- updated_data$ms1_peaks
              ms1_features <- updated_data$ms1_features
              
              }
          }else {
            
            idx <- paste0(ms1_peaks$peakid, "-1") %in% peaks_xcms$xic_label  # all detected peaks
            idx <- which(idx)
            updated_data <- update_df(ms1_peaks, ms1_features, idx, iteration_number = 0)
            ms1_peaks <- updated_data$ms1_peaks
            ms1_features <- updated_data$ms1_features
            
          }
        } else {
          ms1_features[feature_idx, iteration := 0] 
        }
      
      feature_idx <- which(is.na(ms1_features$iteration))[1]
      
      if( is.na(feature_idx) ){
        break  }
    }
    
    # save results
    if( isTRUE(d.out) ){
      dir.create("./Results/")
      
      file_name <- tools::file_path_sans_ext(basename(file_info$mzml_path[sample_idx]))  
      dir.create(paste0("./Results/", file_name))
      saveRDS(ms1_peaks, paste0("./Results/", file_name, '/ms1_peaks.rds'))
      saveRDS(ms1_features, paste0("./Results/", file_name, '/ms1_features.rds'))
      saveRDS(features.l, paste0("./Results/", file_name, '/features_list.rds'))
    }
     
    return(features.l)
  })
  
  names(res) <- tools::file_path_sans_ext(basename(file_info$mzml_path[sample_idx]))
  return(res)
}



#' Update df iteration columns
#'
#' @param ms1_peaks `data.frame` obtains from `DIANMF::extract_xcms_peaks()`
#' @param ms1_features `data.frame` obtains from `DIANMF::extract_xcms_features()`
#' @param idx `c(numerics)` row indexes
#' @param iteration_number `numeric(1)`
#'
#' @return `list(ms1_peaks, ms1_features)`
#' @export
update_df <- function(ms1_peaks, ms1_features, idx, iteration_number){

  # update ms1_peaks
  ms1_peaks$iteration[idx] <- ifelse(
    is.na(ms1_peaks$iteration[idx]),
    iteration_number,
    paste0(ms1_peaks$iteration[idx], ",", iteration_number)
  )
  
  # update ms1_feature
  ms1_features[, iteration := as.character(iteration)]
  ms1_features[sapply(peakidx, function(p) any(unlist(p) %in% idx)), 
               iteration := ifelse(is.na(iteration), as.character(iteration_number), 
                                   paste0(iteration, ",", iteration_number))]
  
  return(list(ms1_peaks = ms1_peaks, ms1_features = ms1_features))
}
