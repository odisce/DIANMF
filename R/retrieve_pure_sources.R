#' Retrieve pure sources: elution profiles and there spectra
#'
#' @param W 
#' @param H 
#' @param ms_type 
#'
#' @return list of pure sources
#' @export 
#' 
#' @import dplyr
pure_sources <- function(W, H, ms_type = c("max", "mean", "sum")){
  res <- lapply(1:ncol(W), function(s){
    # s <- 1
    spect <- as.matrix(W[, s], ncol = 1)
    elut <- as.matrix(H[s, ], nrow = 1) 
    source <- spect %*% t(elut)
    colnames(source) <- colnames(H)
    rownames(source) <- rownames(W)
    
    source <- prepare_mixed_data(ms_mixed = source, mz_values = as.numeric(rownames(source)), rts = as.numeric(colnames(source)) )
    
    if( ms_type == "sum" ){  # Sum of intensities at the same mz
      
      source_spect <- source %>%
        group_by(mz_value) %>%
        summarise(intensity = sum(intensity, na.rm = TRUE), .groups = "drop")
      
    } else if( ms_type == "mean" ){ # mean of intensities at the same mz
      
      source_spect <- source %>%
        group_by(mz_value) %>%
        summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop")
      
    } else {  # max of intensities at the same mz
      
      source_spect <- source %>%
        group_by(mz_value) %>%
        summarise(intensity = max(intensity, na.rm = TRUE), .groups = "drop")
      
    }
    
    return( list( source_eic = source,
                  source_spect = source_spect  ))
  })
}


plot_eics <- function(eics_mat, ms_level = "MS1"){
  
  eics_mat$mz_value <- paste0(ms_level, eics_mat$mz_value)
  
  p <- ggplot2::ggplot( data = eics_mat, aes(x = rt, y = intensity, color = mz_value)) + 
    geom_vline( xintercept = rt_prec, colour = "red", linetype = 2, linewidth = 0.5 ) +
    geom_line() +
    geom_point() +
    xlim( min(eics_mat$rt), max(eics_mat$rt) ) +
    guides(color = 'none')
  
  # if(!is.null(prec_eic_row)){
  #   p <- p +  geom_line(data = eics_mat[eics_mat$mz_value == paste0(ms_level, prec_eic_row), ], aes(x = rt, y = intensity), color = "black")
  # }
  
  return(p)
}


plot_spectrum <- function(pure_spect){
  
  rank <- length(pure_spect)
  W <- matrix(0, ncol = rank, nrow = length(pure_spect[[1]]$mz_value))
  for (i in 1:rank) {
    W[,i] <- pure_spect[[i]]$intensity
  }
  rownames(W) <- pure_spect[[1]]$mz_value
  colnames(W) <- c(paste0("comp",seq(1,rank,1)))
  W <- reshape2::melt(W)
  colnames(W) <- c("mz_value", "comp_nb", "intensity")
  
  p <- ggplot2::ggplot( ) + 
    geom_linerange(data = W, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity, color = comp_nb)) +
    facet_grid(comp_nb~.) +
    theme_bw()
  
  return(p)
}


plot_peak_info <- function(ms_mixed = ms1_mat, W = W_ms1, H = H_ms1, mz_prec, ms_level = c("MS1", "MS2"), pure_sources.l){
  
  mz_ions <- as.numeric(rownames(ms_mixed))
  closest_row <- which.min(abs(mz_ions - mz_prec))
  prec_eic_row <- mz_ions[closest_row]
  
  # plot the mixed elution profiles 
  mixed_data <- prepare_mixed_data(ms_mixed = ms_mixed, mz_values = as.numeric(rownames(ms_mixed)), rts = as.numeric(colnames(ms_mixed)) )
  ms1_mixed_eics <- mixed_data
  ms1_mixed_eics$mz_value <- paste0(ms_level, ms1_mixed_eics$mz_value)
  p_mixed_eics <- plot_eics(eics_mat = ms1_mixed_eics)
  
  # plot mixed spectrum
  ms1_mixed_spectrum <- mixed_data
  p_mixed_spect <- ggplot( ) + 
    geom_linerange(data = ms1_mixed_spectrum, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
    theme_bw()
  
  # plot pure sources
  p_pure_eics <- lapply(1:ncol(W), function(s){
    p_eics <- plot_eics(eics_mat = pure_sources.l[[s]]$source_eic,  ms_level = "MS1")
  })
  p_pure_eics <- patchwork::wrap_plots(p_pure_eics, ncol = 1)  
  
  pure_spect <- lapply(1:ncol(W), function(s){
    res <- pure_sources.l[[s]]$source_spect
  })
  p_pure_spect <- plot_spectrum(pure_spect)
  
  p_mixed <- patchwork::wrap_plots(list(p_mixed_eics, p_mixed_spect), ncol = 2)
  p_pure <-  patchwork::wrap_plots(list(p_pure_eics, p_pure_spect), ncol = 2)
  
  p_all <- patchwork::wrap_plots(list(p_mixed, p_pure), ncol = 1)
  
  return(p_all)
}

