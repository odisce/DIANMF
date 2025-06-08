#' Export MS Spectra
#'
#' @inheritParams get_elutionprofile
#' @param log2L `logic` is `TRUE` for logarithmic transformation, else `FALSE`.
#'
#' @returns `ggplot` of the mixed matrix for the asked feature.
#' @export
#'
#' @import ggplot2 data.table magrittr
#' @importFrom stringr str_to_sentence
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
  targ_feat <- features.l[[feat_coord$sample]]$ms1_features[featureid == feature_id, ][which.max(into), peakid]
  
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
#' @inheritParams plot_EluProfile
#'
#' @returns `ggplot` of the mixed spectra for the asked feature.
#' @export
#'
#' @import ggplot2 data.table magrittr
#' @importFrom stringr str_to_sentence
plot_Spectra <- function(
    features.l = features.l,
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
  targ_feat <- features.l[[feat_coord$sample]]$ms1_features[featureid == feature_id, ][which.max(into), peakid]
  
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

#' Multi-plot mixed/pure feature
#'
#' @inheritParams plot_EluProfile
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


# nev <- get_ms1_rtdiff(msexp_idx) * 1.5
# scan_rt_ext <- 10
# rt_method <- "constant"
# rt_range ?!

#' Plot the XCMS peaks of one iteration after post-deconvolution (i.e. filtering pure sources) 
#'
#' @inheritParams get_spectra
#' @param iteration_index `numeric` iteration index. 
#'
#' @returns ggplot2 plot of peaks kept in the post-deconvolution step.
#' @export 
#' 
#' @import ggplot2
plot_xcms_peaks_range <- function(summary_dt, sample_index, iteration_index){
  
  xcms_peaks <- summary_dt[ sample == sample_index & iteration == iteration_index,  ]
  rt_range <- c( xcms_peaks[1, ]$rtmin-10 , xcms_peaks[1, ]$rtmax+10) 
  
  ggplot(xcms_peaks, aes(rt, mz, group = peakid)) +
    geom_pointrange(aes(x = rt, xmin = rtmin, xmax = rtmax, color = peakfull)) +
    geom_vline(xintercept = rt_range, linetype = 2, color = "red") +
    theme_bw() +
    labs(
      title = "MS1 peak in range",
      caption = sprintf(
        "
            Total MS1 peaks: %i\n
            Peak include in range: %i
          ",
        xcms_peaks[, .N],
        xcms_peaks[peakfull == TRUE, .N]
      )
    )
}

#' Plot rt window with the factorization results
#'
#' @inheritParams get_elutionprofile 
#' @param iteration_idx `numeric` iteration index.
#' @param targets `list` of `data.frame`s of two columns c(compound, mz) of co-eluting compounds, to highlight their eics in the mixed chromatograms plot.
#'
#' @returns ggplot  
#' @export
#' 
#' @import data.table ggplot2 dplyr 
#' @importFrom ggpubr ggarrange
#' @importFrom cowplot ggdraw draw_plot 
plot_rtWind_info <- function(features.l, sample_index = 1, iteration_idx = 1, targets = NULL){
  
  sample_features <- features.l[[sample_index]]$PureFeatures
  # feat_peak_tb <- features.l[[sample_index]]$ms1_features_peaks
  # feat_peak_tb <- feat_peak_tb[sample == sample_index, ]
  info <- sample_features[[iteration_idx]]
  ms_info <- info$ms_info
  
  # plot the pure sources 
  ms1_pure_spectra <- info$ms1_pure_spectra
  temp_x <- ms_info[ms_info$xic_label %in% ms1_pure_spectra$xic_label, ]
  temp_x <- temp_x[, c("xic_label", "mz", "rt")]
  ms1_pure_spectra <- ms1_pure_spectra %>%
    left_join(temp_x, by = "xic_label")
  ms1_pure_spectra <- ms1_pure_spectra[scan_norm >0, ]
  ms1_pure_spectra <- ms1_pure_spectra[value >0, ]
  ms1_pure_spectra$mslevel <- "MS1"
  
  ms2_pure_spectra <- info$ms2_pure_spectra
  temp_x2 <- ms_info[ms_info$xic_label %in% ms2_pure_spectra$xic_label, ]
  temp_x2 <- temp_x2[, c("xic_label", "mz", "rt")]
  ms2_pure_spectra <- ms2_pure_spectra %>%
    left_join(temp_x2, by = "xic_label")
  ms2_pure_spectra <- ms2_pure_spectra[value >0, ]
  ms2_pure_spectra$mslevel <- "MS2"
  
  pure_spectra <- rbind(ms1_pure_spectra, ms2_pure_spectra)
  setDT(pure_spectra)
  pure_spectra[, value := value / max(value), by = .(rank, mslevel)]
  pure_spectra[, value := ifelse(mslevel == "MS1", value, -value)]
  
  p1 <- ggplot(pure_spectra, aes(x = mz, ymin = 0, ymax = value, color = mslevel)) +
    geom_linerange() +
    facet_grid(rank ~ ., scales = "free_y") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    # scale_color_manual(
    #   name = "mslevel",
    #   values = c("MS1" = "black", "MS2" = "darkred"),
    #   labels = c("MS1", "MS2")
    # ) +
    labs(title = "Pure MS1 vs. MS2 spectra", y = "Intensity", x = "m/z") +
    theme_bw() +
    guides(color = 'none')
  
  ms1_pure_chrom <- info$MS1_pure_elution_profiles
  ranks_good <- unique(ms1_pure_chrom$rank) 
  rt_axis <- ms1_pure_chrom[rank == ranks_good[1], ]$rtime
  
  ms1_pure_chrom$mslevel <- "MS1"
  ms1_pure_chrom <- ms1_pure_chrom %>%
    group_by(rank) %>%
    mutate(value = value / max(value, na.rm = TRUE))
  
  # ms2_pure_chrom <- info$MS2_pure_elution_profiles
  # ms2_pure_chrom$mslevel <- "MS2"
  # pure_chrom <- rbind(ms1_pure_chrom, ms2_pure_chrom)
  p2 <- ggplot() +
    geom_line(data = ms1_pure_chrom, aes(x = rtime, y = value), color = "black") +
    # geom_line(data = pure_chrom[pure_chrom$mslevel == "MS2", ], aes(x = rtime, y = -value, color = factor(rank))) +
    facet_grid(rank ~ ., scales = "free_y") +
    labs(
      title = "Pure elution profiles",
      x = "Retention Time (RT)",
      y = "Intensity",
      color = "m/z"
    ) +
    theme_bw() +
    guides(color = 'none')
  
  p_pure <- ggpubr::ggarrange(plotlist = list(p2,p1), ncol = 2)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # plot mixed chromatograms and spectra
  chrom_mix_ms1 <- info$MS1_mixed_mat
  chrom_mix_ms2 <- info$MS2_mixed_mat
  
  ms1_peaks_names <- rownames(chrom_mix_ms1)
  # ms1_peaks_names <- sub("-1$", "", ms1_peaks_names)
  ms1_mz_ions <-  ms_info[ ms_info$xic_label %in% ms1_peaks_names, ]$mz
  # ms1_rt_ions <- ms_info[ ms_info$xic_label %in% ms1_peaks_names, ]$rt
  
  ms2_peaks_names <- rownames(chrom_mix_ms2)
  # ms2_peaks_names <- sub("-2$", "", ms2_peaks_names)
  ms2_mz_ions <-  ms_info[ ms_info$xic_label %in% ms2_peaks_names, ]$mz
  # ms2_rt_ions <- ms_info[ ms_info$xic_label %in% ms2_peaks_names, ]$rt
  
  # mz_ions <- as.numeric(rownames(chrom_mix_ms1))
  # closest_row <- which.min(abs(mz_ions - mz_prec))
  # prec_eic_row <- mz_ions[closest_row]
  
  rownames(chrom_mix_ms1) <- ms1_mz_ions
  colnames(chrom_mix_ms1) <- rt_axis
  chrom_mix_ms1 <- reshape2::melt(chrom_mix_ms1)
  colnames(chrom_mix_ms1) <- c('mz_value', 'rt', 'intensity')
  chrom_mix_ms1$mslevel <- "MS1"
  
  rownames(chrom_mix_ms2) <- ms2_mz_ions
  colnames(chrom_mix_ms2) <- rt_axis
  chrom_mix_ms2 <- reshape2::melt(chrom_mix_ms2)
  colnames(chrom_mix_ms2) <- c('mz_value', 'rt', 'intensity')
  chrom_mix_ms2$mslevel <- "MS2" 
  
  mixed_data <- rbind(chrom_mix_ms1, chrom_mix_ms2)
  
  # Calculate y-axis range for positioning
  y_max1 <- max(mixed_data[mixed_data$mslevel == "MS1", ]$intensity, na.rm = TRUE)
  y_max2 <- -max(mixed_data[mixed_data$mslevel == "MS2", ]$intensity, na.rm = TRUE)
  
  if( !is.null(targets) ){
    closest_matches <- targets %>%
      rowwise() %>%
      mutate(
        closest_mz = {
          ms1_data <- mixed_data %>% filter(mslevel == "MS1")
          ms1_data$mz_value[which.min(abs(ms1_data$mz_value - mz))]
        }
      ) %>%
      ungroup()
    mixed_data$highlight <- ifelse(mixed_data$mz_value %in% closest_matches$closest_mz, TRUE, FALSE)
    annot_data <- mixed_data %>%
      filter(mz_value %in% closest_matches$closest_mz) %>%
      group_by(mz_value) %>%
      slice_max(intensity, n = 1) %>%
      left_join(closest_matches, by = c("mz_value" = "closest_mz"))
    
    p3.1 <- ggplot() +
      # Regular lines (colored by mz)
      geom_line(
        data = mixed_data[!mixed_data$highlight, ],
        aes(x = rt, y = ifelse(mslevel == "MS1", intensity, -intensity), color = factor(mz_value)),
        alpha = 0.6
      ) +
      # Highlighted lines in black
      geom_line(
        data = mixed_data[mixed_data$highlight, ],
        aes(x = rt, y = ifelse(mslevel == "MS1", intensity, -intensity), group = mz_value),
        color = "black", linewidth = 0.5
      ) +
      geom_text(
        data = annot_data,
        aes(x = rt, y = ifelse(mslevel == "MS1", intensity, -intensity), label = compound),
        hjust = -0.1, vjust = 0.5, angle = 90, color = "black", fontface = "plain"
      ) +
      labs(
        title = "Raw chromatograms MS1 vs. MS2",
        x = "Retention Time (RT)",
        y = "Intensity",
        color = "m/z"
      ) +
      theme_bw() +
      guides(color = 'none')
    
    eic_data <- mixed_data %>%
      filter(mz_value %in% closest_matches$closest_mz)
    eic_data <- eic_data %>%
      left_join(closest_matches, by = c("mz_value" = "closest_mz"))
    p3.2 <- ggplot(eic_data, aes(x = rt, y = intensity, color = compound)) +
      geom_line(size = 1) +
      labs(
        # title = "",
        x = "RT",
        y = "Intensity",
        color = "Compound"
      ) +
      theme_minimal()
    
    p3 <- ggdraw() +
      draw_plot(p3.1) +  # Base plot
      draw_plot(p3.2, x = 0.6, y = 0.6, width = 0.38, height = 0.38)  # Inset plot (x, y are from 0 to 1)
  }else{
    p3 <- ggplot() +
      geom_line(
        data = mixed_data,
        aes(x = rt, y = ifelse(mslevel == "MS1", intensity, -intensity), color = factor(mz_value)),
        alpha = 0.6
      ) +
      labs(
        title = "Raw chromatograms MS1 vs. MS2",
        x = "Retention Time (RT)",
        y = "Intensity",
        color = "m/z"
      ) +
      theme_bw() +
      guides(color = 'none')
  }
  
  mixed_data <- mixed_data %>%
    group_by(mz_value, mslevel) %>%
    summarise(intensity = sum(intensity, na.rm = TRUE), .groups = "drop") # ms_type = c("sum", "max", "mean")[1]
  
  # Calculate y-axis range for positioning
  y_max1 <- max(mixed_data[mixed_data$mslevel == "MS1", ]$intensity, na.rm = TRUE)
  y_max2 <- -max(mixed_data[mixed_data$mslevel == "MS2", ]$intensity, na.rm = TRUE)
  
  p4 <- ggplot() +
    geom_linerange(data = mixed_data[mixed_data$mslevel == "MS1", ], aes(x = mz_value, y = intensity, ymin = 0, ymax = intensity, color = factor(mz_value))) +
    geom_linerange(data = mixed_data[mixed_data$mslevel == "MS2", ], aes(x = mz_value, y = -intensity, ymax = 0, ymin = -intensity, color = factor(mz_value))) +
    labs(
      title = "Raw MS1 vs. MS2 spectra",
      x = "m/z",
      y = "Intensity"
    ) +
    annotate("text", x = Inf, y = y_max1 * 0.9, label = "MS1", hjust = 1.1, vjust = 0, fontface = "bold", color = "black") +
    annotate("text", x = Inf, y = y_max2 * 0.9, label = "MS2", hjust = 1.1, vjust = 1, fontface = "bold", color = "black")+ 
    theme_bw() +
    guides(color = 'none')
  
  
  p_mix <- ggpubr::ggarrange(plotlist = list(p3,p4), ncol = 2)
  
  p_final <- ggpubr::ggarrange(plotlist = list(p_mix,p_pure), 
                               nrow = 2, 
                               heights = c(1, 2))
  
  # ggsave(paste0(d.plot, '/ff.png'), p_final, w=30, h = 15, dpi = 300)
  
  pure_spectra[mslevel == "MS2", value := abs(value)]
  
  return(list(
    "pure_spectra" = pure_spectra,
    "p_final" = p_final
  ))
} 
