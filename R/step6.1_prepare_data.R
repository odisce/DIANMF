#' Prepare mixed MS1 and MS2 data.
#'
#' @param ms_mixed mixed `matrix` of a peak. Every row is an extracted ion chromatogram (EIC) along the peak rt scans.
#' @param mz_values Fragments mz values of the eics in ms_mixed. The row TRUE names.
#' @param rts EICs retention time values in ms_mixed. The columns TRUE names.
#'
#' @return `data.frame` for every EIC with the mz and rt information.
#' 
#' @export
#' 
#' @importFrom reshape2 melt
prepare_mixed_data <- function(ms_mixed = mat1, mz_values = as.numeric(rownames(ms1_mat)), rts = as.numeric(colnames(ms1_mat)) ){
  
  rownames(ms_mixed) <- mz_values
  colnames(ms_mixed) <- rts
  ms_mixed <- reshape2::melt(ms_mixed)
  colnames(ms_mixed) <- c('mz_value', 'rt', 'intensity')
  
  return(ms_mixed)
}

#-------------------------------------------------------------------------------

#' Prepare pure spectra for plotting.
#'
#' @param W Pure spectra `matrix`.
#'
#' @return Pure components spectra information `data.frame`.
#' 
#' @export
#' 
#' @importFrom reshape2 melt
#' @import magrittr
prepare_pure_spectra <- function(W){
  nb_components <- ncol(W)
  colnames(W) <- c(paste0("comp",seq(1,nb_components,1)))
  W <- reshape2::melt(W)
  colnames(W) <- c("mz_value", "comp_nb", "intensity")
  W <- W %>% 
    group_by(comp_nb) %>%
    mutate(intensity = intensity / max(intensity)) 
  
  return(W)
}

#-------------------------------------------------------------------------------

#' Prepare pure elution profiles (chromatograms) for plotting.
#'
#' @param H Elution profiles `matrix`.
#' @param rt `numeric` vector chromatograms retention time values in H.
#'
#' @return Pure components chromatograms information `data.frame`.
#' 
#' @export
#' 
#' @importFrom reshape2 melt
prepare_pure_eics <- function(H, rt){
 
  rownames(H) <- paste0("comp", seq(1,nrow(H)))
  colnames(H) <- rt
  H <- reshape2::melt(H)
  colnames(H) <- c("comp_nb", "rt", "intensity")
  
  return(H)
}

#-------------------------------------------------------------------------------

#' Plot mixed and pure elution profiles.
#'
#' @inheritParams prepare_mixed_data
#' @param ms_pure_H Pure elution profiles `data.frame`. 
#' @inheritParams extract_mixedMat
#' @inheritParams check_ms1_ions
#' @inheritParams choose_ms2_pure_spectrum
#'
#' @return ggplot2 plot
#' 
#' @export
#' 
#' @import ggplot2
#' @import patchwork
plot_MS_eics <- function(ms_mixed, ms_pure_H = NULL, ms_level = c("MS1", "MS2"), rt_prec, choosen_comp){
  
  # plot the MIXED elution profiles 
  ms1_mixed_eics <- prepare_mixed_data(ms_mixed = ms_mixed, mz_values = as.numeric(rownames(ms_mixed)), rts = as.numeric(colnames(ms_mixed)) )
  ms1_mixed_eics$mz_value <- paste0(ms_level, ms1_mixed_eics$mz_value)
  p1 <- ggplot2::ggplot( data = ms1_mixed_eics, aes(x = rt, y = intensity, color = mz_value)) + 
    geom_vline( xintercept = rt_prec, colour = "red", linetype = 1, linewidth = 0.5 ) +
    geom_line() +
    geom_point() +
    xlim( min(ms1_mixed_eics$rt), max(ms1_mixed_eics$rt) ) +
    guides(color = 'none')
  
  
  if( is.null(ms_pure_H) ){
    return(p1)
  } else {
    # plot the PURE elution profiles
    ms_rt <- as.numeric(colnames(ms_mixed))
    ms1_pure_eics <- prepare_pure_eics(H = ms_pure_H, rt = ms_rt)
    p2 <- ggplot2::ggplot(data = ms1_pure_eics, aes(rt, intensity , color = comp_nb)) +
      facet_grid(comp_nb~.) +
      geom_line() +
      geom_point() +
      xlim(min(ms1_pure_eics$rt), max(ms1_pure_eics$rt)) +
      labs( caption = paste(ms_level, " elution profiles; rt:", round(rt_prec,4), " good comp:", choosen_comp))
    
    p <- p1 / p2 + plot_layout(ncol = 1, guides = "collect")  # Align plots vertically with shared legend
    
    return(p)
  }
 
}

#-------------------------------------------------------------------------------

#' plot mixed and pure spectra
#'
#' @inheritParams prepare_mixed_data
#' @param ms_pure_W `data.frame` of pure spectra
#' @inheritParams plot_MS_eics
#' @inheritParams filter_ms2_spectrum
#'
#' @return ggplot2 plot
#' 
#' @export
#' 
#' @import ggplot2
plot_MS_spectra <- function(ms_mixed, ms_pure_W, ms_level = c("MS1", "MS2"), mz_prec, choosen_comp){
  
  # plot the MIXED MS1 spectrum
  ms1_mixed_spectrum <- prepare_mixed_data(ms_mixed = ms_mixed, mz_values = as.numeric(rownames(ms_mixed)), rts = as.numeric(colnames(ms_mixed)) )
  ms1_mixed_spectrum$intensity <- ms1_mixed_spectrum$intensity / max(ms1_mixed_spectrum$intensity)
  ms1_mixed_spectrum$comp_nb <- "comp0-mixed"
  ms1_mixed_spectrum <- ms1_mixed_spectrum[, c(1,3,4)]
  
  # plot PURE MS1 spectra
  ms1_pure_spectra <- prepare_pure_spectra(W = ms_pure_W)
  
  spectra <- rbind(ms1_mixed_spectrum, ms1_pure_spectra)
  
  p <- ggplot2::ggplot( ) + 
    geom_linerange(data = spectra, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity, color = comp_nb)) +
    facet_grid(comp_nb~.) +
    labs( caption = paste(ms_level, " spectra; mz:", round(mz_prec,4), " good_comp:", choosen_comp)) +
    theme_bw();
  
  return(p)
}

#-------------------------------------------------------------------------------

#' Plot measured vs library spectra.
#'
#' @param measured_spectrum `data.frame` measured spectrum.
#' @param library_spectrum `data.frame` library/reference spectrum.
#'
#' @return ggplot2 plot
#' 
#' @export
#' 
#' @import ggplot2
plot_spectra_vs <- function(measured_spectrum, library_spectrum){
  measured_spectrum <- as.data.frame(measured_spectrum)
  measured_spectrum$intensity <- measured_spectrum$intensity / max(measured_spectrum$intensity)
  
  library_spectrum <- as.data.frame(library_spectrum)
  colnames(library_spectrum) <- c('mz_value', 'intensity')
  library_spectrum$intensity <- library_spectrum$intensity / max(library_spectrum$intensity)
  
  p <- ggplot2::ggplot( ) +
    geom_linerange(data = measured_spectrum, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
    geom_linerange(data = library_spectrum, stat = "identity", aes( x = mz_value, y = -intensity, ymin = -intensity, ymax = 0, color = 'red')) +
    guides(color = FALSE) +
    labs( caption = "measured VS library spectra.") +
    theme_bw();
  
  return(p)
}

#-------------------------------------------------------------------------------

# plot_all_info <- function(d.plot,
#                           peak.idx,
#                           ms_mixed1, ms_mixed2,
#                           ms_pure_H1, ms_pure_H2,
#                           ms_pure_W1, ms_pure_W2,
#                           rt_prec, mz_prec,
#                           choosen_comp_ms1 , choosen_comp_ms2){
# 
#   if (!dir.exists(d.plot)) {
#     dir.create(d.plot, recursive = TRUE)
#   };
#   
#   # for MS2 data
#   ms1_elution_profiles <- plot_MS_eics(ms_mixed1, ms_pure_H1, ms_level = "MS1", rt_prec = rt_prec, choosen_comp = choosen_comp_ms1)
#   ms1_spectra <- plot_MS_spectra(ms_mixed1, ms_pure_W1, ms_level = "MS1", mz_prec, choosen_comp = choosen_comp_ms1)
#   
#   # for MS2 data
#   ms2_elution_profiles <- plot_MS_eics(ms_mixed2, ms_pure_H2, ms_level = "MS2", rt_prec = rt_prec, choosen_comp = choosen_comp_ms2)
#   ms2_spectra <- plot_MS_spectra(ms_mixed2, ms_pure_W2, ms_level = "MS2", mz_prec, choosen_comp = choosen_comp_ms2)
#   
#   ggsave(ms1_elution_profiles, filename = paste0(d.plot, '/', peak.idx, '_ms1_elution_profiles.svg'), device = "svg");
#   ggsave(ms1_spectra, filename = paste0(d.plot, '/', peak.idx, '_ms1_spectra.svg'), device = "svg");
#   ggsave(ms2_elution_profiles, filename = paste0(d.plot, '/', peak.idx, '_ms2_elution_profiles.svg'), device = "svg");
#   ggsave(ms2_spectra, filename = paste0(d.plot, '/', peak.idx, '_ms2_spectra.svg'), device = "svg");
# 
# }

#-------------------------------------------------------------------------------