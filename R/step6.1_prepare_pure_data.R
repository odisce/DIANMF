# This script contains functions to reshape\plot pure and mixed EICS and spectra.

# I want to make many changes in these functions.!!!!!!!!!!!!!!!!

#' prepare mixed eics for plotting
#'
#' @param eics description
#' @param retention_times description
#'
#' @importFrom reshape2 melt
#'
#' @return data.frame
#'
prepare_mixed_eics <- function(eics, retention_times){

  eics <- melt(eics)
  colnames(eics) <- c("MS_type", "scan_idx_nb", "intensity")

  v <- c()
  v <- lapply(1:nrow(eics), function(i){

    #i <- 1
    v[i] <- retention_times[eics$scan_idx[i]]

  })
  v <- as.numeric(v)
  eics <- cbind(eics, rt=v)
  # eics$intensity <- eics$intensity / max(eics$intensity)

  return(eics)
}

#-------------------------------------------------------------------------------

#' prepare pure eics for plotting
#'
#' @param H description
#' @param retention_times description
#'
#' @importFrom reshape2 melt
#'
#' @return data.frame
#'
prepare_pure_eics <- function(H, retention_times){
  nb_components <- nrow(H)
  H <- t(H)
  colnames(H) <- c(paste0("comp",seq(1,nb_components,1)))
  H <- melt(H)
  colnames(H) <- c("scan_idx", "comp_nb", "intensity")

  rt_2 <- c()
  rt_2 <- lapply(1:nrow(H), function(i){

    rt_2[i] <- retention_times[H$scan_idx[i]]

  })
  rt_2 <- as.numeric(rt_2)
  class(rt_2)
  H <- cbind(H, rt_2)
  colnames(H) <- c("scan_idx", "comp_nb", "intensity", "rt")
  #H$intensity <- H$intensity / max(H$intensity)

  return(H)
}

#-------------------------------------------------------------------------------

#' prepare pure spectra for plotting
#'
#' @param W description
#' @param mz_values description
#'
#' @importFrom reshape2 melt
#'
#' @return data.frame
#'
prepare_pure_spectra<- function(W , mz_values){
  nb_components <- ncol(W)
  colnames(W) <- c(paste0("comp",seq(1,nb_components,1)))
  W <- melt(W)
  colnames(W) <- c("MS_scan", "comp_nb", "intensity")

  v <- c()
  for (i in 1:nrow(W)) {

    v[i] <- mz_values[[ W$MS_scan[i] ]]
  }
  W <- cbind( W, v )

  colnames(W) <- c("MS_scan", "comp_nb", "intensity", "mz_value")
  #W$intensity <- W$intensity / max(W$intensity)

  return(W)
}

#-------------------------------------------------------------------------------

#' prepare mixed spectra for plotting
#'
#' @param spectra description
#' @param mz_values description
#'
#' @importFrom reshape2 melt
#'
#' @return data.frame
#'
prepare_mixed_spectra <- function(spectra, mz_values){

  spectra <- melt(spectra)
  colnames(spectra) <- c("MS_type", "scan_idx_nb", "intensity")

  v <- c()
  for (i in 1:nrow(spectra)) {

    v[i] <- mz_values[ spectra$scan_idx_nb[i] ]
  }

  spectra <- cbind(spectra, v)
  colnames(spectra)[4] <- 'mz_values'

  df_max <- spectra %>%
    group_by(MS_type) %>%
    summarize(max_value = max(intensity))

  df_max <- cbind(df_max, mz_values)

  colnames(df_max) <- c('MS_type', 'intensity', 'mz_value')

  return(df_max)
}

#-------------------------------------------------------------------------------

#' plot the mixed eics, experimental and real rt of the precursor
#'
#' @param data.f1 data.frame
#' @param rt.n_prec numeric precursor apex real rt (the rt I have in the compounds table)
#' @param rt.n_peak numeric peak apex rt (peak apex rt detected by xcms)
#'
#' @importFrom ggplot2 ggplot
#'
#' @return ggplot2 plot of mixed eics
#'
plot_mixed_eics <- function(data.f1, rt.n_prec, rt.n_peak){

  plot1 <- ggplot( data = data.f1, aes(x = rt, y = intensity, color = MS_type)) +
    geom_vline( xintercept = as.numeric( rt.n_prec ), colour = "black", linetype = 1, size = 0.5 ) +
    geom_vline( xintercept = as.numeric( rt.n_peak ), colour = "red", linetype = 2, size = 0.5 ) +
    geom_point(alpha = 0.5) +
    geom_line() +
    geom_line(data = subset(data.f1, MS_type == "MS1"), color = "black", size = 0.8, linetype = 2) +
    xlim( min(data.f1$rt), max(data.f1$rt) ) +
    guides(color = FALSE)  +
    theme(
      panel.background = element_rect(fill = "white"),  # Sets the background color to white
      panel.grid.major = element_blank(),  # Removes major gridlines
      panel.grid.minor = element_blank(),  # Removes minor gridlines
      panel.border = element_blank(),  # Removes panel border
      axis.line = element_line(color = "black") ) +
    ggtitle( paste('prec_rt =', round(rt.n_prec,2), 'ms1_peak_rt=', round(rt.n_peak,2), 'rt_diff= ', round(abs(rt.n_prec - rt.n_peak),2) ))
}

#-------------------------------------------------------------------------------
