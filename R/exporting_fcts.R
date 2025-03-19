# the following function, extract for every identified pure source in all iterations one Spectra object;
# this object contains 1 MS1 pure spectra and the 10 MS2 pure spectra from different isolation windows.

#' Export MS Spectra
#'
#' @param features.l `list` obtainer with `DIANMF::DIANMF.f()`
#' @param max_method `string` to choose a method to select the optimal iteration:
#'   - `contribution`: select the first iteration with maximum contribution.
#'   - `max_value`: select the first iteration with maximum apex pure value.
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table magrittr
get_feature_summary <- function(
  features.l,
  max_method = c("contribution", "max_value")[1]
) {
  feature_dt <- data.table()
  for (spli in seq_len(length(features.l))) {
    for (iteri in seq_len(length(features.l[[spli]]$PureFeatures))) {
      if (is.null(features.l[[spli]]$PureFeatures[[iteri]])) {
        next
      }
      out_dt <- features.l[[spli]]$PureFeatures[[iteri]]$ms_info[!is.na(source),]
      out_dt[, iteration := as.integer(iteri)]
      feature_dt <- rbind(feature_dt, out_dt)
      rm(out_dt)
    }
  }
  if (max_method == "max_value") {
    output <- feature_dt[order(-apex_val), .SD[which.max(apex_val), ], by = .(featureid, sample)]
  } else if (max_method == "contribution") {
    output <- feature_dt[order(-apex_val), .SD[which.max(contribution), ], by = .(featureid, sample)]
  } else {
    stop(sprintf("max_method not recognized: %s", max_method))
  }
  
  return(output[])
}


#' Export MS Spectra
#'
#' @inheritParams get_elutionprofile
#'
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table
get_feature_coord <- function(
  features.l,
  summary_dt = NULL,
  feature_id,
  sample_index = NULL,
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l, max_method = max_method)
  }
  feat_coord <- summary_dt[featureid == feature_id, ]
  if (!is.null(sample_index)) {
    feat_coord <- feat_coord[sample == sample_index, ]
  }
  if (nrow(feat_coord) == 0) {
    message(sprintf("Feature %s not extracted in sample: %i", feature_id, sample_index))
    return(FALSE)
  }
  return(feat_coord)
}

#' Extract elution profile
#'
#' @param type `string` to retrieve the "pure" or "mixed" elution profiles.
#' @param method `string` to extract "all" sources or just the "best".
#' @inheritParams get_feature_coord
#'
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table magrittr
get_elutionprofile <- function(
  features.l = features,
  summary_dt = NULL,
  feature_id,
  sample_index = NULL,
  type = c("pure", "mixed")[1],
  method = c("all", "best")[2],
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l, max_method = max_method)
  }
  feat_coord <- get_feature_coord(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    max_method = max_method
  )
  if (isFALSE(feat_coord)) {
    return(FALSE)
  }

  ## Extract MS1 / MS2
  if (type == "mixed") {
    mixed_names <- c("MS1_mixed_mat", "MS2_mixed_mat")
    rt_names <- c("MS1_pure_elution_profiles", "MS2_pure_elution_profiles")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i], rt_names[i]) %in% i_names)) {
        mixed_mat_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]] %>%
          melt() %>%
          as.data.table()
        rt_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[rt_names[i]]][, .(scan_norm, rtime)] %>%
          unique()
        setnames(mixed_mat_dt, c("xic_label", "scan_norm", "value"))
        output <- merge(mixed_mat_dt, rt_dt, by = "scan_norm")
        output[, mslevel := mslve_vc[i]]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
    mixed_mat_out <- mixed_mat_out[, .(xic_label, rtime, value, mslevel, rank = as.integer(NA))]
  } else if (type == "pure") {
    mixed_names <- c("MS1_pure_elution_profiles", "MS2_pure_elution_profiles")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i]) %in% i_names)) {
        source_to_get <- features[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][, unique(rank)]
        if (method == "best") {
          source_to_get <- feat_coord[, source]
        }
        output <- features[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][rank %in% source_to_get,]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
    mixed_mat_out <- mixed_mat_out[, .("xic_label" = rank, rtime, value, "mslevel" = MSid, rank = rank)]
  } else {
    stop('type argument not recognized')
  }
  return(mixed_mat_out)
}

#' Extract Spectra
#'
#' @inheritParams get_feature_coord
#' @inheritParams get_elutionprofile
#'
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table magrittr
get_spectra <- function(
  features.l,
  summary_dt = NULL,
  feature_id,
  sample_index = 1,
  type = c("pure", "mixed")[1],
  method = c("all", "best")[2],
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l, max_method = max_method)
  }
  feat_coord <- get_feature_coord(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    max_method = max_method
  )
  if (isFALSE(feat_coord)) {
    return(FALSE)
  }
  ## Extract MS1 / MS2
  if (type == "mixed") {
    mixed_names <- c("MS1_mixed_mat", "MS2_mixed_mat")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i], "ms_info") %in% i_names)) {
        mixed_mat_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]] %>%
          melt() %>%
          as.data.table()
        setnames(mixed_mat_dt, c("xic_label", "scan_norm", "value"))
        mixed_mat_dt <- mixed_mat_dt[, .(value = max(value)), by = .(xic_label)]
        ms_info_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][["ms_info"]]
        output <- merge(
          mixed_mat_dt,
          ms_info_dt[, .(xic_label, msLevel, IsoWin, featureid, mz, rt, rank = as.integer(NA))],
          by = "xic_label"
        )
        output[, mslevel := mslve_vc[i]]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
  } else if (type == "pure") {
    mixed_names <- c("ms1_pure_spectra", "ms2_pure_spectra")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i]) %in% i_names)) {
        source_to_get <- features[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][, unique(rank)]
        if (method == "best") {
          source_to_get <- feat_coord[, source]
        }
        output <- features[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][rank %in% source_to_get,]
        ms_info_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][["ms_info"]]
        output <- merge(
          output,
          ms_info_dt[, .(xic_label, msLevel, IsoWin, featureid, mz, rt)],
          by = "xic_label"
        )
        output[, mslevel := mslve_vc[i]]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
  } else {
    stop('type argument not recognized')
  }
  
  return(mixed_mat_out[, .(xic_label, IsoWin, mslevel, mz, value, rank)])
}

#' Export MS Spectra
#'
#' @param features.l `list` obtainer with `DIANMF::DIANMF.f()`
#' @param summary_dt
#' @param feature_id
#' @param sample_index
#' @param log2L
#' @inheritParams 
#'
#' @returns `ggplot` of the mixed matrix for the asked feature.
#' @export
#'
#' @import ggplot2 data.table magrittr
plot_EluProfile <- function(
  features.l,
  summary_dt = NULL,
  feature_id = NULL,
  sample_index = 1,
  log2L = FALSE,
  type = c("pure", "mixed")[1],
  method = c("all", "best")[2],
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l)
  }
  feat_coord <- summary_dt[featureid == feature_id & sample == sample_index,]
  mixed_mat_out <- get_elutionprofile(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    type = type,
    max_method = max_method,
    method = method
  )

  if (log2L) {
    mixed_mat_out[, x := log2(value + 1), by = .(xic_label)]
  } else {
    mixed_mat_out[, x := value]
  }
  if (mixed_mat_out[, length(unique(mslevel))] == 2) {
    mixed_mat_out[, x := ifelse(grepl("MS1", mslevel), x, -x)]
  }

  ## Get feature related peaks
  targ_feat <- features[[feat_coord$sample]]$ms1_features[featureid == feature_id, ][which.max(into), peakid]

  ggplot(mixed_mat_out, aes(rtime, x, group = interaction(xic_label, mslevel), color = mslevel)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.6) +
    geom_line(data = mixed_mat_out[
      rtime %between% feat_coord[, c(rtmin, rtmax)]
    ][grepl(targ_feat, xic_label), ], alpha = 1, color = "red") +
    theme_bw() +
    facet_grid(rank ~ .) +
    scale_y_continuous(labels = function(x) sprintf("%0.2g", x)) +
    labs(
      title = sprintf(
        "%s matrix of feature: %s",
        stringr::str_to_sentence(type),
        feature_id
      ),
      x = "Retention time (in sec.)",
      y = "Intensity"
    )
}

#' Plot MS Spectra
#'
#' @param features.l `list` obtainer with `DIANMF::DIANMF.f()`
#' @param summary_dt
#' @param feature_id
#' @param sample_index
#' @param log2L
#'
#' @returns `ggplot` of the mixed spectra for the asked feature.
#' @export
#'
#' @import ggplot2 data.table magrittr
plot_Spectra <- function(
  features.l = features,
  summary_dt = NULL,
  feature_id = NULL,
  sample_index = 1,
  log2L = FALSE,
  type = c("pure", "mixed")[1],
  method = c("all", "best")[2],
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l)
  }
  feat_coord <- summary_dt[featureid == feature_id & sample == sample_index,]
  mixed_mat_out <- get_spectra(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    type = type,
    max_method = max_method,
    method = method
  )

  if (log2L) {
    mixed_mat_out[, x := log2(value + 1), by = .(xic_label)]
  } else {
    mixed_mat_out[, x := value]
  }
  if (mixed_mat_out[, length(unique(mslevel))] == 2) {
    mixed_mat_out[, x := ifelse(grepl("MS1", mslevel), x, -x)]
  }

  ## Get feature related peaks
  targ_feat <- features[[feat_coord$sample]]$ms1_features[featureid == feature_id, ][which.max(into), peakid]

  ggplot(mixed_mat_out, aes(
    x = mz,
    y = x,
    ymin = 0,
    ymax = x,
    group = interaction(xic_label, mslevel),
    color = as.factor(IsoWin)
  )) +
    geom_hline(yintercept = 0) +
    geom_linerange(alpha = 0.6) +
    geom_linerange(data = mixed_mat_out[
        grepl(targ_feat, xic_label),
      ],
      alpha = 1,
      color = "red"
    ) +
    theme_bw() +
    facet_grid(rank ~ .) +
    scale_y_continuous(labels = function(x) sprintf("%0.2g", x)) +
    labs(
      title = sprintf(
        "%s matrix of feature: %s",
        stringr::str_to_sentence(type),
        feature_id
      ),
      x = "m/Z",
      y = "Intensity"
    )
}

#' Multiplot mixed/pure feature
#'
#' @inheritParams plot_Spectra
#'
#' @returns `ggplot` of the mixed spectra for the asked feature.
#' @export
#'
#' @import ggplot2 data.table magrittr
plot_feature <- function(
  features.l = features,
  summary_dt = NULL,
  feature_id = NULL,
  sample_index = 1,
  log2L = FALSE,
  max_method,
  method = c("all", "best")[2]
) {
  feat_coord <- get_feature_coord(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    max_method = max_method
  )
  if (isFALSE(feat_coord)) {
    warning(
      sprintf(
        "features %s not successfully extracted in sample: %s",
        feature_id,
        sample_index
      )
    )
    return(NULL)
  }
  plotA <- plot_EluProfile(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    type = "mixed",
    log2L = log2L,
    max_method = max_method,
    method = method
  )
  plotB <- plot_EluProfile(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    type = "pure",
    log2L = log2L,
    max_method = max_method,
    method = method
  )
  plotC <- plot_Spectra(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    type = "mixed",
    log2L = log2L,
    max_method = max_method,
    method = method
  )
  plotD <- plot_Spectra(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    type = "pure",
    log2L = log2L,
    max_method = max_method,
    method = method
  )
  ggpubr::ggarrange(
    ggpubr::ggarrange(plotA, plotC, ncol = 2, align = "hv", common.legend = TRUE),
    ggpubr::ggarrange(plotB, plotD, ncol = 2, align = "hv", common.legend = TRUE),
    nrow = 2,
    align = 'hv',
    common.legend = TRUE
  )
}
#' Export MS Spectra
#'
#' @param features.l `list`
#'
#' @returns `list` of `Spectra` objects
#' @export
#' 
#' @importFrom Spectra Spectra
#' @import dplyr
exportMSSpectra <- function(features.l) {
  
  MSSpectra <- lapply(seq_along(features.l), function(iteration) {
    
    ms1PureSpectra <- features.l[[iteration]]$ms1_pure_spectra %>%
      arrange(rank) %>%
      group_by(rank) %>%
      group_split()
    
    ms2PureSpectra <- features.l[[iteration]]$ms2_pure_spectra %>%
      arrange(isolationWindowTargetMz) %>%
      group_by(rank) %>%
      group_split()
    
    info <- features.l[[iteration]]$xcms_assigned_sources %>%
      arrange(source) %>%
      group_by(source) %>%
      group_split()
    
    res <- lapply(seq_along(ms1PureSpectra), function(source_idx) {
      # MS1
      ms1_df <- ms1PureSpectra[[source_idx]]
      ms1_spd <- S4Vectors::DataFrame(
        msLevel = 1L,
        name = iteration,
        id = unique(ms1_df$rank),
        rtime = mean(info[[source_idx]]$rt),
        COMMENT = paste(info[[source_idx]]$xic_label, collapse = ";")
      )
      
      ms1_spd$mz <- list(ms1_df$mz)
      ms1_spd$intensity <- list(ms1_df$value)
      ms1_spectra <- Spectra::Spectra(ms1_spd)
      
      # MS2
      ms2_df <- ms2PureSpectra[[source_idx]]
      isolation_spectra <- lapply(unique(ms2_df$isolationWindowTargetMz), function(targetMz) {
        subset_df <- ms2_df[ms2_df$isolationWindowTargetMz == targetMz, ]  
        
        ms2_spd <- S4Vectors::DataFrame(
          msLevel = 2L,
          name = iteration,
          id = unique(subset_df$rank), 
          isolationWindowTargetMz = as.numeric(targetMz),  
          isolationWindowLowerMz = as.numeric(unique(subset_df$isolationWindowLowerMz)),
          isolationWindowUpperMz = as.numeric(unique(subset_df$isolationWindowUpperMz))  )
        
        ms2_spd$mz <- list(subset_df$mz)
        ms2_spd$intensity <- list(subset_df$value)
        return(Spectra::Spectra(ms2_spd))
      })
      
      isolation_spectra <- do.call(c, isolation_spectra)
      msSpectra <- c(ms1_spectra, isolation_spectra)
      return(msSpectra) 
    })
    
    return(res)
  })
  
  MSSpectra <- unlist(MSSpectra, recursive = FALSE)
  return(MSSpectra)
}