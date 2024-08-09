#' prepare mixed MS1 and MS2 data
#'
#' @param m matrix
#' @param mz_value numeric vector contains the fragments mz values of the eics in m
#' @param rt numeric vector contains the fragments retention time values of the eics in m
#'
#' @return data.frame 
#' @export
#' @importFrom reshape2 melt
prepare_mixed_data <- function(m = mat1, mz_value = as.numeric(rownames(ms1_mat)), rt = as.numeric(colnames(ms1_mat)) ){
  
  rownames(m) <- mz_value
  colnames(m) <- rt
  m <- reshape2::melt(m)
  colnames(m) <- c('mz_value', 'rt', 'intensity')
  
  return(m)
}

#-------------------------------------------------------------------------------

#' prepare pure spectra for plotting
#'
#' @param W matrix contains pure spectra
#'
#' @return data.frame
#' @export
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

#' prepare pure elution profiles for plotting
#'
#' @param H matrix contains elution profiles
#' @param rt numeric vector contains the fragments retention time values of the eics in H
#'
#' @return data.frame
#' @export
#' @importFrom reshape2 melt
prepare_pure_eics <- function(H, rt){
 
  rownames(H) <- paste0("comp", seq(1,nrow(H)))
  colnames(H) <- rt
  H <- reshape2::melt(H)
  colnames(H) <- c("comp_nb", "rt", "intensity")
  
  return(H)
}

#-------------------------------------------------------------------------------

#' Select the good MS1 spectrum related to the MS1 peak
#'
#' @param W_ms1 matrix of MS1 pure spectra
#' @param mz_prec numeric(1) mz value of the precursor
#'
#' @return data.frame pure MS1 spectrum
#' @export
choose_ms1_pure_spectrum <- function(W_ms1, mz_prec){
  
  ms1_pure_spectra <- prepare_pure_spectra(W_ms1);
  
  mz_ms1_ions <- round(as.numeric(rownames(W_ms1)), 4);
  closest_row <- which.min(abs(mz_ms1_ions - mz_prec));
  choosen_comp_ms1 <- which.max(W_ms1[closest_row, ]);
  
  ms1_pure_spectrum <- ms1_pure_spectra[ms1_pure_spectra$comp_nb == paste0('comp',choosen_comp_ms1), ];
  ms1_pure_spectrum <- ms1_pure_spectrum[, c('mz_value', 'intensity')];
  ms1_pure_spectrum <- ms1_pure_spectrum[ms1_pure_spectrum['intensity'] != 0, ];
  
  return(list(
    'ms1_pure_spectrum' = ms1_pure_spectrum,
    'choosen_comp_ms1' = choosen_comp_ms1  ))
}

#-------------------------------------------------------------------------------

#' Select the good MS2 spectrum related to the MS1 peak
#'
#' @param W_ms2 matrix of MS2 pure spectra
#' @param choosen_comp_ms2 numeric(1) good component
#'
#' @return data.frame pure MS2 spectrum
#' @export
#' @import magrittr
choose_ms2_pure_spectrum <- function(W_ms2, choosen_comp_ms2){
  
  ms2_pure_spectra <- prepare_pure_spectra(W_ms2);
  
  ms2_pure_spectra <-  ms2_pure_spectra %>%
    group_by(comp_nb) %>%
    mutate(intensity = intensity / max(intensity));
  
  ms2_pure_spectrum <- subset(ms2_pure_spectra, comp_nb == paste0('comp',choosen_comp_ms2));
  ms2_pure_spectrum <- ms2_pure_spectrum[, c('mz_value','intensity')];
  
  return(ms2_pure_spectrum)
}

#-------------------------------------------------------------------------------

#' plot mixed and pure elution profiles
#'
#' @param ms_pure_H data.frame of pure elution profiles 
#' @param ms_level character MS1 or MS2
#' @param rt_prec numeric precursor retention time
#' @param ms_mixed matrix
#' @param choosen_comp integer the index of the good component
#'
#' @return ggplot2 plot
#' @export
#' @import ggplot2
#' @import ggpubr
plot_MS_eics <- function(ms_mixed, ms_pure_H, ms_level = c("MS1", "MS2"), rt_prec, choosen_comp){
  
  # plot the MIXED elution profiles 
  ms1_mixed_eics <- prepare_mixed_data(m = ms_mixed, mz_value = as.numeric(rownames(ms_mixed)), rt = as.numeric(colnames(ms_mixed)) )
  ms1_mixed_eics$mz_value <- paste0(ms_level, ms1_mixed_eics$mz_value)
  p1 <- ggplot2::ggplot( data = ms1_mixed_eics, aes(x = rt, y = intensity, color = mz_value)) + 
    geom_vline( xintercept = rt_prec, colour = "red", linetype = 1, size = 0.5 ) +
    geom_line() +
    geom_point() +
    xlim( min(ms1_mixed_eics$rt), max(ms1_mixed_eics$rt) ) +
    guides(color = FALSE)
  
  # plot the PURE elution profiles
  ms_rt <- as.numeric(colnames(ms_mixed))
  ms1_pure_eics <- prepare_pure_eics(H = ms_pure_H, rt = ms_rt)
  p2 <- ggplot2::ggplot(data = ms1_pure_eics, aes(rt, intensity , color = comp_nb)) +
    facet_grid(comp_nb~.) +
    geom_line() +
    geom_point() +
    xlim(min(ms1_pure_eics$rt), max(ms1_pure_eics$rt)) +
    labs( caption = paste(ms_level, " elution profiles; rt:", round(rt_prec,4), " good comp:", choosen_comp))
  
  p <- ggpubr::ggarrange( p1, p2, ncol =  1, nrow = 2, align = "h", common.legend = TRUE )
  
  return(p)
}

#-------------------------------------------------------------------------------

#' plot mixed and pure spectra
#'
#' @param ms_mixed matrix
#' @param ms_pure_W data.frame of pure spectra
#' @param ms_level character MS1 or MS2
#' @param mz_prec numeric precursor mz
#' @param choosen_comp integer the index of the good component
#'
#' @return ggplot2 plot
#' @export
#' @import ggplot2
plot_MS_spectra <- function(ms_mixed, ms_pure_W, ms_level = c("MS1", "MS2"), mz_prec, choosen_comp){
  
  # plot the MIXED MS1 spectrum
  ms1_mixed_spectrum <- prepare_mixed_data(m = ms_mixed, mz_value = as.numeric(rownames(ms_mixed)), rt = as.numeric(colnames(ms_mixed)) )
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

#' return 4 ggplots of MS1 and MS2 mixed and pure data of a specific peak
#'
#' @param d.plot path to save the plots
#' @param ms_mixed1 ms1 mixed matrix
#' @param ms_mixed2 ms2 mixed matrix
#' @param ms_pure_H1 ms1 pure elution profiles matrix
#' @param ms_pure_H2 ms2 elution profiles matrix
#' @param ms_pure_W1 ms1 pure spectra matrix
#' @param ms_pure_W2 ms2 pure spectra matrix
#' @param rt_prec prec/peak rt
#' @param mz_prec prec/peak mz
#' @param choosen_comp_ms1 the good ms1 comp
#' @param choosen_comp_ms2 the good ms2 comp
#' @param peak.idx integer peak index
#'
#' @return 4 ggplot2 plots that contains all MS1, MS2 mixed and pure information
plot_all_info <- function(d.plot,
                          peak.idx,
                          ms_mixed1, ms_mixed2,
                          ms_pure_H1, ms_pure_H2,
                          ms_pure_W1, ms_pure_W2,
                          rt_prec, mz_prec,
                          choosen_comp_ms1 , choosen_comp_ms2){

  if (!dir.exists(d.plot)) {
    dir.create(d.plot, recursive = TRUE)
  };
  
  # for MS2 data
  ms1_elution_profiles <- plot_MS_eics(ms_mixed1, ms_pure_H1, ms_level = "MS1", rt_prec = rt_prec, choosen_comp = choosen_comp_ms1)
  ms1_spectra <- plot_MS_spectra(ms_mixed1, ms_pure_W1, ms_level = "MS1", mz_prec, choosen_comp = choosen_comp_ms1)
  
  # for MS2 data
  ms2_elution_profiles <- plot_MS_eics(ms_mixed2, ms_pure_H2, ms_level = "MS2", rt_prec = rt_prec, choosen_comp = choosen_comp_ms2)
  ms2_spectra <- plot_MS_spectra(ms_mixed2, ms_pure_W2, ms_level = "MS2", mz_prec, choosen_comp = choosen_comp_ms2)
  
  ggsave(ms1_elution_profiles, filename = paste0(d.plot, '/', peak.idx, '_ms1_elution_profiles.svg'), device = "svg");
  ggsave(ms1_spectra, filename = paste0(d.plot, '/', peak.idx, '_ms1_spectra.svg'), device = "svg");
  ggsave(ms2_elution_profiles, filename = paste0(d.plot, '/', peak.idx, '_ms2_elution_profiles.svg'), device = "svg");
  ggsave(ms2_spectra, filename = paste0(d.plot, '/', peak.idx, '_ms2_spectra.svg'), device = "svg");

}

#-------------------------------------------------------------------------------
