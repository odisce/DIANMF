#' Process xcms object
#'
#' @inheritParams extract_xcms_peaks
#' @param dir_out `logical` TRUE to save results, else FALSE
#' @param sample_idx NULL to process all files else, `numeric` index of the sample.
#' @param MS2_ISOEACHL `logical`
#' @param MS1MS2_L `logical`
#' @inheritParams get_spectra_values
#' @inheritParams nGMCAs
#' @param rt_method c("peak", "constant")
#' @inheritParams has_n_consecutive_non_zero
#' @param scan_rt_ext `numeric`
#' @param min_distance `numeric`
#' @param featuresn `numeric` to limit the number if iteration (k) on each file for testing purpose, set it to NULL to extract all features.
#' @param verbose `logical` to show execution messages
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
DIANMF.f <- function(
  msexp,
  dir_out = FALSE,
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
  rt_method = c("peak", "constant"),
  scan_rt_ext = 10,
  min_distance = 5,
  nscans = 4,
  featuresn = NULL,
  combineSpectra_arg = list(peaks = "intersect", ppm = 4, tolerance = 0.005, minProp = 0.05),
  verbose = FALSE
) {
  options(datatable.verbose = FALSE)
  ms1_peaks <- extract_xcms_peaks(msexp)
  # ms1_peaks[, iteration := as.character(NA)]
  ms1_features_all <- extract_xcms_features(msexp, quantifyL = TRUE)
  ms1_features_peaks <- ms1_features_all[, .(peakindex = unlist(peakidx)), by = featureid]
  ms1_features_peaks <- merge(ms1_peaks, ms1_features_peaks, by = "peakindex")
  file_info <- MsExperiment::sampleData(msexp) %>% as.data.table()

  # Check and Filter the file_info based on MSExperiment indexes
  if (!is.null(sample_idx)) {
    if (any(sample_idx > length(msexp))) {
      stop("sample_idx doesn't match msexp object indexes.")
    } else {
      file_info <- file_info[sample_idx,]
    }
  }

  ## Iterate over files ----
  res <- lapply(seq_len(file_info[, .N]), function(s_idx) {
    # s_idx = 1
    ms1_features <- copy(ms1_features_peaks[sample == s_idx, ])
    ms1_features[, iteration := as.character(NA) ]
    msexp_idx <- xcms::filterFile(msexp, s_idx)
    nev <- get_ms1_rtdiff(msexp_idx) * 1.5
    s_idx_name <- file_info[s_idx, basename(mzml_path)]
    if (verbose) {
      message(
        sprintf(
          "Processing sample: %i %s",
          s_idx,
          s_idx_name
        )
      )
    }
    min_rt <- ms1_features[, min(rtmin)]
    max_rt <- ms1_features[, max(rtmax)]
    features.l <- list()
    feature_idx <- 1
    k <- 0
    ## Sort and Extract features
    ms1_features <- ms1_features[order(-into), ]
    while (feature_idx <= nrow(ms1_features)) {
      # Option to limit the number of features to extrac
      k <- k + 1
      if (!is.null(featuresn)) {
        if (k > featuresn) {
          break
        }
      }
      # Set current feature to process
      feature_idx <- ms1_features[is_filled == 0 & is.na(iteration), ][1, featureid]
      if (is.na(feature_idx)) {
        break
      }
      if (verbose) {
        message(sprintf("feature index: %i ------ k: %i", feature_idx, k))
      }
      peak_i <- ms1_features[featureid == feature_idx & is_filled == 0, ][which.max(into), ]
      if (nrow(peak_i) <= 0) {
        if (verbose) {
          message("    No peak found for this feature: skipping")
        }
        ms1_features[featureid == feature_idx, iteration := 0]
        next
      } else {
        rt_range <- peak_i[which.max(into), range(c(rtmin, rtmax)) + c(-scan_rt_ext - nev, +scan_rt_ext + nev)]
        if (rt_method == "constant") {
          rt_range <- (mean(rt_range) + c(-scan_rt_ext - nev, +scan_rt_ext + nev))
        }
        rt_range[1] <- max(min_rt, rt_range[1])  # to avoid rt < 0
        rt_range[2] <- min(max_rt, rt_range[2])  # to avoid rt out of range
        ## Subset msexp object
        msexp_idx_rt <- xcms::filterRt(msexp_idx, rt_range)
        ## Generate MS1 peaklist
        if (verbose) {
          message("    Generating MS1 peaks list")
        }
        ms1_peaks_i <- ms1_features[rtmin <= rt_range[2] & rtmax >= rt_range[1], ]
        ## flag peaks partially, fully or apex inside the window
        ms1_peaks_i[, peakfull := ifelse(
          (rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
          ifelse(rt %between% rt_range, "apex", "partial")
        )]
        ms1_peaks_i[, msLevel := 1]
        ms1_peaks_i[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
        ## Align scans and get time dictionary
        if (verbose) {
          message("    Aligning MS(n) events")
        }
        time_dic <- align_scans(
          msexp = msexp_idx_rt,
          rt_range = rt_range,
          sample_idx = 1
        )
        ## Get MS2 peaklist
        if (verbose) {
          message("    Generating MS2 peaks list")
        }
        ms2_peaks_i <- generate_peaklist(
          msexp = msexp_idx_rt,
          mslevel = 2,
          ms2isowinL = TRUE,
          combineSpectra_arg = combineSpectra_arg,
          ppm = 6
        )
        ## Extract raw data
        if (verbose) {
          message("    Extracting ions signals")
        }
        raw_dt_i <- get_spectra_values(msexp_idx_rt)
        ## Combine MS1 MS2 peaklists if asked
        if (MS1MS2_L) {
          peaks_ls <- list(
            "MS1MS2" = rbind(
              ms1_peaks_i,
              ms2_peaks_i,
              fill = TRUE
            )
          )
        } else {
          peaks_ls <- list(
            "MS1" = ms1_peaks_i,
            "MS2" = ms2_peaks_i
          )
        }
        ## NMF
        ms2H <- NULL ## If MS2 after MS1, pass H to init MS2
        nmf_result_ls <- list(
          "ms_info" = NULL,
          "mixed_mat" = NULL,
          "pure_elution_profiles" = NULL,
          "pure_spectra" = NULL
        )
        NextIter <- FALSE
        for (i in seq_len(length(peaks_ls))) {
          if (verbose) {
            message(sprintf("    Generating %s XICs", peaks_ls[[i]][, .N]))
          }
          ms_xics_i <- build_XICs(
            peaks_dt = peaks_ls[[i]],
            rawdt = raw_dt_i
          )
          ms_xics_i <- merge(ms_xics_i, time_dic[, .(rtime, scan_norm)], by = "rtime")
          if (verbose) {
            message("    Build mixed matrix")
          }
          ms_mixed_i <- build_mixed_matrix(ms_xics_i, nscans = nscans)
          if (nrow(ms_mixed_i$mixedmat) <= 0) {
            if (verbose) {
              message("    Empty mixed matrix: skipping")
            }
            NextIter <- TRUE
            break
          }
          if (verbose) {
            message(sprintf("    Unmixing %s", names(peaks_ls)[i]), appendLF = FALSE)
            timeA <- Sys.time()
          }
          rank <- min(rank, (ncol(ms_mixed_i$mixedmat) - 1))
          nmf_i <- nGMCAs(
            X.m = ms_mixed_i$mixedmat,
            rank = rank,
            maximumIteration = maximumIteration,
            maxFBIteration = maxFBIteration,
            toleranceFB = toleranceFB,
            initialization_method = initialization_method,
            H_sub = ms2H,
            errors_print = errors_print,
            method = method
          )
          if (verbose) {
            timeB <- Sys.time()
            message(sprintf(": Done (%0.2fs)", difftime(timeB, timeA, units = "secs")), appendLF = TRUE)
            message("    Extracting pure sources", appendLF = FALSE)
          }
          # Format output
          pure_rt <- reshape2::melt(nmf_i$A) %>% as.data.table()
          setnames(pure_rt, c("rank", "scan_norm", "value"))
          pure_rt[, MSid := names(peaks_ls)[i]]
          ## Add mean rtime
          pure_rt <- merge(
            pure_rt,
            time_dic[, .(rtime = mean(rtime)), by = .(scan_norm)],
            by = "scan_norm"
          )
          pure_mz <- reshape2::melt(nmf_i$S) %>% as.data.table()
          setnames(pure_mz, c("xic_label", "rank", "value"))
          pure_mz[, MSid := names(peaks_ls)[i]]
          ## Add contribution
          pure_mz[, contribution := value / sum(value), by = xic_label]
          ## Add max source contribution to peaks
          xic_max_contrib <- pure_mz[, .(source = .SD[which.max(contribution), rank]), by = .(xic_label)]
          peaks_ls[[i]] <- merge(
            peaks_ls[[i]],
            xic_max_contrib,
            by = "xic_label"
          )
          if (verbose) {
            timeC <- Sys.time()
            message(sprintf(": Done (%0.2fs)", difftime(timeC, timeB, units = "secs")), appendLF = TRUE)
          }
          nmf_result_ls <- list(
            "ms_info" = rbind(
              nmf_result_ls$ms_info,
              peaks_ls[[i]],
              fill = TRUE
            ),
            "mixed_mat" = rbind(
              nmf_result_ls$mixed_mat,
              ms_mixed_i$mixedmat
            ),
            "pure_rt" = rbind(
              nmf_result_ls$pure_rt,
              pure_rt
            ),
            "pure_spectra" = rbind(
              nmf_result_ls$pure_spectra,
              pure_mz
            )
          )
          if (i == 1 && "MS1" %in% names(peaks_ls)[i] && length(peaks_ls) > 1) {
            # Passing H to next iteration
            if (verbose) {
              message(sprintf("    Passing H to %s", names(peaks_ls)[i + 1]))
            }
            ms2H <- nmf_i$A
          } else {
            ms2H <- NULL
          }
          rm(ms_mixed_i, pure_rt, pure_mz, ms_xics_i, nmf_i)
        }
        if (NextIter) {
          ms1_features[featureid == feature_idx, iteration := ifelse(is.na(iteration), 0, iteration)]
          next
        }
        ## Filter sources based on MS1 contribution threshold
        tempdt <- nmf_result_ls$pure_spectra[xic_label %in% nmf_result_ls$ms_info[msLevel == 1, unique(xic_label)],]
        good_sources <- tempdt[, .(contL = max(contribution) >= 0.6), by = .(rank, MSid)][contL == TRUE, .(rank, MSid)]
        if (good_sources[, .N] <= 0) {
          if (verbose) {
            message("    no source passing threshold, skipping")
          }
          ms1_features[featureid == feature_idx, iteration := ifelse(is.na(iteration), 0, iteration)]
          next
        }
        ## Detect pure elution profile peaks
        pure_rt_good <- nmf_result_ls$pure_rt[good_sources, on = c("rank", "MSid")]
        chromM <- list()
        for (ranki in pure_rt_good[, seq_len(length(unique(rank)))]) {
          chromM[[ranki]] <- MSnbase::Chromatogram(
            pure_rt_good[rank == unique(rank)[ranki]][order(rtime), rtime],
            pure_rt_good[rank == unique(rank)[ranki]][order(rtime), value]
          )
        }
        pure_peaks <- pure_rt_good[order(rtime), {
          output <- MSnbase::Chromatogram(rtime, value) %>% {
              MSnbase::MChromatograms(list(.))
            } %>%
            xcms::findChromPeaks(
              object = .,
              param =  CentWaveParam(
                ppm = 0,
                prefilter = c(3, 1),
                peakwidth = c(2, 10),
                snthresh = 1,
                mzCenterFun = "wMeanApex3",
                integrate = 1,
                mzdiff = 0,
                noise = 0,
                firstBaselineCheck = FALSE)
            ) %>%
            xcms::chromPeaks() %>%
            as.data.table()
          if (output[, .N] <= 0) {
            NULL
          } else {
            output
          }
        }, by = .(rank, MSid)]

        if (pure_peaks[, .N] <= 0) {
          if (verbose) {
            message("    No peaks detected in any sources: skipping")
          }
          ms1_features[featureid == feature_idx, iteration := ifelse(is.na(iteration), 0, iteration)]
          next
        }
        ## Filter detected peaks and select the two most intense
        pure_peaks_f <- pure_peaks[sn >= 0, lapply(.SD[order(-into),], head, 2), by = .(rank, MSid)]
        ## Remove peaks to close to the border unless at the file rt limits
        peak_lims <- rt_range
        peak_lims[1] <- ifelse(peak_lims[1] <= min_rt, peak_lims[1] - 3, peak_lims[1] + 3)
        peak_lims[2] <- ifelse(peak_lims[2] >= max_rt, peak_lims[2] + 3, peak_lims[2] - 3)
        if (peak_lims[1] <= min_rt) {peak_lims[1] <- peak_lims[1] - 10} else {}
        if (peak_lims[2] >= max_rt) {peak_lims[2] <- peak_lims[2] + 10} else {}
        pure_peaks_f <- pure_peaks_f[rt %between% peak_lims, ]
        if (pure_peaks_f[, .N] <= 0) {
          if (verbose) {
            message("    No peaks remains in any sources: skipping")
          }
          ms1_features[featureid == feature_idx, iteration := ifelse(is.na(iteration), 0, iteration)]
          next
        }
        ## Attribute MS1 peaks to each good sources
        for (PPeak_row in pure_peaks_f[, seq_len(.N)]) {
          pk_range <- pure_peaks_f[PPeak_row, c(rtmin, rtmax)]
          rank_sel <- pure_peaks_f[PPeak_row, rank]
          nmf_result_ls$ms_info[
            msLevel == 1L &
              source == rank_sel &
              rt %between% pk_range & peakfull %in% c("apex", "full"),
            iteration := k
          ]
        }
        features.l[[k]] <- list(
          "MS1_mixed_mat" = nmf_result_ls$mixed_mat[
            rownames(nmf_result_ls$mixed_mat) %in% nmf_result_ls$ms_info[msLevel == 1, xic_label],
          ],
          "MS2_mixed_mat" = nmf_result_ls$mixed_mat[
            rownames(nmf_result_ls$mixed_mat) %in% nmf_result_ls$ms_info[msLevel == 2, xic_label],
          ],
          "ms1_pure_spectra" = nmf_result_ls$pure_spectra[
            xic_label %in% nmf_result_ls$ms_info[msLevel == 1, xic_label],
          ],
          "ms2_pure_spectra" = nmf_result_ls$pure_spectra[
            xic_label %in% nmf_result_ls$ms_info[msLevel == 2, xic_label],
          ],
          "MS1_pure_elution_profiles" = nmf_result_ls$pure_rt[grepl("MS1", MSid),],
          "MS2_pure_elution_profiles" = nmf_result_ls$pure_rt[grepl("MS2", MSid),],
          "xcms_assigned_sources" = pure_peaks_f,
          "ms2_info" = nmf_result_ls$ms_info[msLevel == 2, ]
        )
        ## Retrieve all features containing extracted peaks
        extracted_peakindex <- ms1_peaks[peakid %in% nmf_result_ls$ms_info[!is.na(iteration), unique(peakid)], unique(peakindex)]
        ## Remove features where peaks are fully included in the current window
        nextracted_peakindex_full <- nmf_result_ls$ms_info[peakfull %in% "full" & !peakindex %in% extracted_peakindex, peakindex]
        extracted_peakindex_all <- c(extracted_peakindex, nextracted_peakindex_full) %>% unique()
        ## even if not attributed (no point in extracting them again)
        ms1_features_peaks[peakindex %in% extracted_peakindex_all,
          iteration := {
            iteration %>%
              strsplit(., ",") %>%
              unlist() %>%
              c(., k) %>%
              unique() %>%
              na.omit() %>%
              paste0(., collapse = ",")
          }
        ]
        feature_extracted <- ms1_features_peaks[peakindex %in% extracted_peakindex_all, unique(featureid)]
        feature_nextracted <- ms1_features_peaks[!peakindex %in% extracted_peakindex & peakindex %in% nextracted_peakindex_full, unique(featureid)]
        ms1_features[
          featureid %in% feature_extracted,
          iteration := {
            iteration %>%
              strsplit(., ",") %>%
              unlist() %>%
              c(., k) %>%
              unique() %>%
              na.omit() %>%
              paste0(., collapse = ",")
          },
          by = .(featureid)
        ]
        if (verbose) {
          message(
            sprintf(
              "    Extracted %i (%i not in sources) features in iteration %i (%i/%i/%i cumulated/remaining/total)",
              length(feature_extracted),
              length(feature_nextracted),
              k,
              ms1_features[!is.na(iteration), .N],
              ms1_features[is.na(iteration), .N],
              ms1_features[, .N]
            )
          )
        }
        ## If targeted feature not in one component, return 0
        ms1_features[featureid == feature_idx, iteration := ifelse(is.na(iteration), 0, iteration)]
      }
    }
    
    # # save results
    # if( isTRUE(dir_out) ){
    #   dir.create("./Results/")
      
    #   file_name <- tools::file_path_sans_ext(basename(file_info$mzml_path[sample_idx]))  
    #   dir.create(paste0("./Results/", file_name))
    #   saveRDS(ms1_peaks, paste0("./Results/", file_name, '/ms1_peaks.rds'))
    #   saveRDS(ms1_features, paste0("./Results/", file_name, '/ms1_features.rds'))
    #   saveRDS(features.l, paste0("./Results/", file_name, '/features_list.rds'))
    # }
     
    return(
      list(
        "PureFeatures" = features.l,
        "ms1_features_peaks" = ms1_features_peaks,
        "ms1_features" = ms1_features
      )
    )
  })
  names(res) <- tools::file_path_sans_ext(basename(file_info$mzml_path[sample_idx]))
  return(res)
}


#' Update df iteration columns
#'
#' @param ms1_peaks `data.frame` obtains from `DIANMF::extract_xcms_peaks()`
#' @param ms1_features `data.frame` obtains from `DIANMF::extract_xcms_features()`
#' @param idx `c(numeric)` row indexes
#' @param iteration_number `numeric(1)`
#'
#' @return A `list` containing updated `ms1_peaks` and `ms1_features` data.tables.
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
