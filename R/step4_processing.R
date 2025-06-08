#' Process xcms object
#' 
#' @inheritParams align_scans
#' @param dir_out `logical` TRUE to save results, else FALSE
#' @param MS2_ISOEACHL `logical`
#' @param MS1MS2_L `logical`
#' @inheritParams get_spectra_values
#' @inheritParams nGMCAs
#' @param rt_method c("peak", "constant")
#' @inheritParams has_n_consecutive_non_zero
#' @param scan_rt_ext `numeric` retention time extension of the rt window.
#' @param min_contrib `numeric` contribution proportion to keep ions (between 0 and 1).
#' @param min_distance `numeric` mark the allowed distance from the rt window border for a peak to be entirely in the window or partially.
#' @param featuresn `numeric` to limit the number if iteration (k) on each file for testing purpose, set it to `NULL` to extract all features.
#' @param verbose `logical` to show execution messages.
#' @param clean_sources `TRUE` to clean/remove sources with less than 0.6 contribution for any peaks, `FALSE` to keep all sources.
#' @param BPPARAM A `BiocParallelParam` object defining the parallelization backend to be used. Defaults to the value returned by `BiocParallel::bpparam()`.
#'
#' @return `list` of info of deconvoluted (pure) sources with the peaks-features data.
#' @export
#' 
#' @importFrom xcms findChromPeaks chromPeaks
#' @importFrom MSnbase MChromatograms Chromatogram
#' @importFrom reshape2 melt
#' @importFrom data.table setnames setDT
#' @import magrittr
#' @import dplyr
#' @import BiocParallel
#' @importFrom MsExperiment sampleData
#' @importFrom tools file_path_sans_ext
DIANMF.f <- function(
  msexp,
  dir_out = FALSE,
  sample_idx = 1,
  MS2_ISOEACHL = TRUE,
  MS1MS2_L = TRUE,
  rank = 10,
  maximumIteration = 200,
  maxFBIteration = 100,
  toleranceFB = 1e-05,
  min_contrib = 0.6,
  initialization_method = "nndsvd",
  errors_print = FALSE,
  method = "svds",
  sparsityA = FALSE,
  rt_method = c("peak", "constant"),
  scan_rt_ext = 10,
  min_distance = 5,
  nscans = 4,
  featuresn = NULL,
  clean_sources = TRUE,
  combineSpectra_arg = list(peaks = "intersect", ppm = 4, tolerance = 0.005, minProp = 0.05),
  verbose = FALSE,
  BPPARAM = bpparam()
) {
  # msexp = xcms_obj
  # dir_out = FALSE
  # sample_idx = 1
  # MS2_ISOEACHL = TRUE
  # MS1MS2_L = TRUE
  # rank = 20
  # min_contrib = 0.6
  # maximumIteration = 200
  # maxFBIteration = 100
  # toleranceFB = 1e-05
  # initialization_method = "nndsvd"
  # errors_print = FALSE
  # method = "svds"
  # scan_rt_ext = 10
  # min_distance = 4
  # featuresn = 2
  # nscans = 6
  # rt_method = "constant"
  # clean_sources = TRUE
  # combineSpectra_arg = list(peaks = "union", ppm = 5, tolerance = 0.005)
  # verbose = TRUE
  # options(datatable.verbose = FALSE)
  ms1_peaks <- extract_xcms_peaks(msexp)
  ms1_features_all <- extract_xcms_features(msexp, quantifyL = TRUE)
  ms1_features_peaks <- ms1_features_all[, .(peakindex = unlist(peakidx), mzmed, rtmed), by = featureid]
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
  res <- BiocParallel::bplapply(
    X = file_info$InjectionOrder,
    FUN = function(s_idx) {
      # s_idx = 1
      ms1_features <- copy(ms1_features_peaks[sample == s_idx, ])
      ms1_features[, iteration := as.character(NA) ]
      msexp_idx <- xcms::filterFile(msexp, s_idx)
      nev <- get_ms1_rtdiff(msexp_idx) * 1.5
      s_idx_name <- file_info[InjectionOrder == s_idx, ]$file_name
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
      last_feature <- ms1_features[, last(featureid)]
      ## Extract features ----
      while (feature_idx != last_feature) {
        # Option to limit the number of features to extrac
        k <- k + 1
        if (!is.null(featuresn)) {
          if (k > featuresn) {
            break
          }
        }
        # Set current feature to process
        feature_idx <- ms1_features[is_filled == 0 & is.na(iteration), ][1, featureid]
        # feature_idx <- "FT04898"
        
        if (is.na(feature_idx)) {
          break
        }
        if (verbose) {
          message(sprintf("feature index: %s ------ k: %i", feature_idx, k))
        }
        peak_i <- ms1_features[featureid == feature_idx & is_filled == 0, ][which.max(into), ]
        if (nrow(peak_i) <= 0) {
          if (verbose) {message("    No peak found for this feature: skipping")}
          next
        } else {
          if (rt_method == "constant") {
            rt_range <- peak_i[, rt + c(-scan_rt_ext - nev, +scan_rt_ext + nev)]
          } else if (rt_method == "peak") {
            rt_range <- peak_i[, range(c(rtmin, rtmax)) + c(-scan_rt_ext - nev, +scan_rt_ext + nev)]
          } else {
            stop(sprintf("rt_method not recognized: %s", rt_method))
          }
          rt_range[1] <- max(min_rt, rt_range[1])  # to avoid rt < 0
          rt_range[2] <- min(max_rt, rt_range[2])  # to avoid rt out of range
          ## Subset msexp object
          msexp_idx_rt <- xcms::filterRt(msexp_idx, rt_range) %>% suppressMessages()
          ## Generate MS1 peaklist
          if (verbose) {
            message("    Generating MS1 peaks list")
          }
          ms1_peaks_i <- ms1_features[rtmin <= rt_range[2] & rtmax >= rt_range[1], -c("iteration")]
          ## flag peaks partially, fully, apex inside apex_border from the window
          border_lim <- min(c(min_distance, floor(diff(rt_range) / 3)))
          rt_limits <- rt_range + c(+border_lim, -border_lim)
          ms1_peaks_i[, peakfull := ifelse((rtmin >= rt_range[1] & rtmax <= rt_range[2]), "full",
            ifelse(!rt %between% rt_limits, "apex_border",
              ifelse(rt %between% rt_range, "apex", "partial"
          )))]
          ms1_peaks_i[, msLevel := 1]
          ms1_peaks_i[, xic_label := paste0(peakid, "-", msLevel), by = .(peakid, msLevel)]
          ms1_peaks_i[, intensity := into]
          features_iter <- ms1_peaks_i[peakfull %in% c("full", "apex"), unique(featureid)]
          ## Exclude features when rtmed not in range
          feat_notinR <- ms1_peaks_i[!rtmed %between% rt_limits, unique(featureid)]
          features_iter <- features_iter %>% {.[!. %in% feat_notinR]} %>% c(., feature_idx) %>% unique()
          ## Set all those features to '0' if not already extracted
          ms1_features[featureid %in% features_iter & is.na(iteration), iteration := "0"]
  
          if (verbose) {
            message(
              sprintf(
                "    Unmixing %i/%i features (peaks: %i/%i/%i, apex+full/border/partial)",
                length(features_iter),
                ms1_features[is.na(iteration), length(unique(featureid))],
                ms1_peaks_i[peakfull %in% c("full", "apex"), length(unique(peakid))],
                ms1_peaks_i[peakfull == "apex_border", length(unique(peakid))],
                ms1_peaks_i[peakfull == "partial", length(unique(peakid))]
              )
            )
          }
          ## Align scans and get time dictionary
          if (verbose) {
            message("    Aligning MS(n) events")
          }
          time_dic <- align_scans(
            msexp = msexp_idx_rt,
            rt_range = rt_range,
            sample_idx = 1
          )
          ## Get normalized time dic
          time_dic_i <- time_dic[, .(rtime = mean(rtime)), by = .(scan_norm)]
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
          ## NMF unmixing ----
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
              timeX <- Sys.time()
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
            rank <- min(rank, (ncol(ms_mixed_i$mixedmat)))
            nmf_i <- nGMCAs(
              X.m = ms_mixed_i$mixedmat,
              rank = rank,
              maximumIteration = maximumIteration,
              maxFBIteration = maxFBIteration,
              toleranceFB = toleranceFB,
              initialization_method = initialization_method,
              H_sub = ms2H,
              errors_print = errors_print,
              method = method,
              sparsityA = sparsityA
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
              time_dic_i,
              by = "scan_norm"
            )
            pure_mz <- reshape2::melt(nmf_i$S) %>% as.data.table()
            setnames(pure_mz, c("xic_label", "rank", "value"))
            pure_mz[, MSid := names(peaks_ls)[i]]
  
            ## Add XCMS peaks contribution to MS1
            if (grepl("MS1", names(peaks_ls)[i])) {
              ### Get scan_norm at apex
              peaks_ls[[i]][peakfull %in% c("full", "apex"), "scan_norm" := {
                time_dic_i[which.min(abs(rtime - rt)), scan_norm]
              }, by = xic_label]
              pure_mz <- merge(
                pure_mz,
                peaks_ls[[i]][, .(xic_label, scan_norm)],
                by = "xic_label"
              )
              pure_mz[!is.na(scan_norm), apex_val := {
                val_vc <- nmf_i$A[rank, scan_norm] * nmf_i$S[xic_label, rank]
                val_vc
              }, by = .(xic_label, rank)]
  
              pure_mz[, contribution := apex_val / sum(apex_val, na.rm = TRUE), by = .(xic_label)]
              pure_mz[is.na(contribution), contribution := 0]
              # pure_mz[, contribution := value / sum(value), by = xic_label]
              ## Add max source contribution to peaks
              xic_max_contrib <- pure_mz[, .SD[which.max(contribution),], by = .(xic_label)][, .(xic_label, source = rank, contribution, apex_val, iteration = k)]
              peaks_ls[[i]] <- merge(
                peaks_ls[[i]],
                xic_max_contrib,
                by = "xic_label",
                all = TRUE
              )
            }
            if (verbose) {
              timeC <- Sys.time()
              message(sprintf(": Done (%0.2fs)", difftime(timeC, timeB, units = "secs")), appendLF = TRUE)
            }
            nmf_result_ls <- list(
              "ms_info" = data.table::rbindlist(
                list(
                  nmf_result_ls$ms_info,
                  peaks_ls[[i]]
                ),
                fill = TRUE
              ),
              "mixed_mat" = rbind(
                nmf_result_ls$mixed_mat,
                ms_mixed_i$mixedmat
              ),
              "pure_rt" = data.table::rbindlist(
                list(
                  nmf_result_ls$pure_rt,
                  pure_rt
                ),
                fill = TRUE
              ),
              "pure_spectra" = data.table::rbindlist(
                list(
                  nmf_result_ls$pure_spectra,
                  pure_mz
                ),
                fill = TRUE
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
            next
          }
          ## Reorder MS2 sources using similarity with MS1
          if (all(c("MS1", "MS2") %in% nmf_result_ls$pure_rt[, unique(MSid)])) {
            ## get ms1 pure rt
            ms1prt <- nmf_result_ls$pure_rt[MSid == "MS1", ] %>% {
                dcast(.[order(scan_norm),], rank ~ scan_norm, value.var = "value")
              } %>% {
                as.matrix(.[, -1])
              }
            ms2prt <- nmf_result_ls$pure_rt[MSid == "MS2", ] %>% {
                dcast(.[order(scan_norm),], rank ~ scan_norm, value.var = "value")
              } %>% {
                as.matrix(.[, -1])
              }
            cor_dt <- data.table()
            for (i in seq_len(nrow(ms1prt))) {
              for (y in seq_len(nrow(ms2prt))) {
                cor_val <- cor(ms1prt[i, ], ms2prt[y, ], method = "pearson")
                outdt <- data.table("rank_1" = i, "rank_2" = y, cor_sc = cor_val)
                cor_dt <- rbind(cor_dt, outdt)
              }
            }
            cor_dt <- cor_dt[order(-cor_sc), ]
            rank_dic <- data.table()
            for (rkn in seq_len(nrow(ms2prt))) {
              outdt <- cor_dt[1, ]
              rank_dic <- rbind(rank_dic, outdt)
              cor_dt <- cor_dt[
                !rank_1 %in% outdt$rank_1 &
                  !rank_2 %in% outdt$rank_2,
              ]
            }
            ## reorder MS1 sources
            new_order <- rank_dic[order(rank_2), rank_1]
            nmf_result_ls$pure_rt[MSid == "MS2", rank := new_order[rank]]
            nmf_result_ls$pure_spectra[MSid == "MS2", rank := new_order[rank]]
          }
          ## Clean sources with less than 0.6 contribution for any peaks
          if (clean_sources) {
            source_to_keep <- nmf_result_ls$ms_info[contribution >= min_contrib, unique(source)]
            nmf_result_ls$pure_rt <- nmf_result_ls$pure_rt[rank %in% source_to_keep,]
            nmf_result_ls$pure_spectra <- nmf_result_ls$pure_spectra[rank %in% source_to_keep,]
            rm(source_to_keep)
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
            "ms_info" = nmf_result_ls$ms_info
          )
  
          ## Mark valid features
          valid_features <- nmf_result_ls$ms_info[featureid %in% features_iter & !is.na(source) & contribution >= min_contrib, unique(featureid)]
          ms1_features[featureid %in% valid_features, iteration := {
            iteration %>%
              strsplit(., ",") %>%
              unlist() %>%
              c(., k) %>%
              unique() %>%
              {.[. != "0"]} %>%
              na.omit(.) %>%
              paste0(., collapse = ",")
          }]
          ## Print message
          if (verbose) {
            timeZ <- Sys.time()
            message(
              sprintf(
                "    iter %i : Extracted %i/%i elligible features %i/%i/%i (cumulated/remaining/total) in: %0.2fs",
                k,
                length(valid_features),
                length(features_iter),
                ms1_features[!is.na(iteration), length(unique(featureid))],
                ms1_features[is.na(iteration), length(unique(featureid))],
                ms1_features[, length(unique(featureid))],
                difftime(timeZ, timeX, units = "secs")
              )
            )
          }
        }
      }
      
      l_output <-  list(
        "PureFeatures" = features.l,
        "ms1_features_peaks" = ms1_features_peaks,
        "ms1_features" = ms1_features
      )
        
      rm(ms1_features)
      return(l_output)
    },
    BPPARAM = BPPARAM
  )
  names(res) <- file_info$file_name
  return(res)
}
