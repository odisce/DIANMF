#' Calculate the elution profiles correlation to choose the good component.
#'
#' @param chromo_main `numeric` MS1 chromatogram.
#' @param chromos `numeric` MS2 chromatograms.
#'
#' @return `numeric(1)` index of the MS2 chromatogram best correlates with the MS1 chromatogram.
#' 
#' @export
#' 
#' @examples
#'  m <- matrix(c(5,0.3,0.5,1,0.9,1,1.5,1.4,1.5,2,1.8,2), nrow = 3, ncol = 4)
#'  x <- matrix(c(0.5,1,1.5,2), nrow = 4, ncol = 1)
#'  
#'  elutions_corelation(chromo_main = x, chromos = m)
#' 
#' @importFrom stats cor
elutions_corelation <- function(chromo_main, chromos){
  correlations <- c()
  for (i in 1:nrow(chromos)) {
    correlations[i] <- stats::cor(chromo_main, chromos[i, ])
  }
  comp <- which.max(correlations)
  return(comp)
}


#' Select the good MS1 spectrum related to the precursor.
#'
#' @param W_ms1 MS1 pure spectra `matrix`.
#' @inheritParams filter_ms2_spectrum
#'
#' @return Precursor MS1 pure spectrum (normalized) `data.frame`.
#'  
#' @export
#' 
#' @examples
#'  W <- matrix(c(0.2, 5.0, 0.6, 0.4, 0.1, 7.0, 0.3, 0.5, 0.3, 6.0, 0.4, 0.2), nrow = 4, byrow = TRUE)
#'  rownames(W) <- c(1,2,3,4)
#'  mz <- 3
#'  
#'  extract_ms1_pure_spectrum(W_ms1 = W, mz_prec = mz)
#' 
extract_ms1_pure_spectrum <- function(W_ms1, mz_prec){
  
  ms1_pure_spectra <- prepare_pure_spectra(W_ms1)
  
  mz_ms1_ions <- round(as.numeric(rownames(W_ms1)), 4)
  closest_row <- which.min(abs(mz_ms1_ions - mz_prec))
  choosen_comp_ms1 <- which.max(W_ms1[closest_row, ])
  
  ms1_pure_spectrum <- ms1_pure_spectra[ms1_pure_spectra$comp_nb == paste0('comp',choosen_comp_ms1), ]
  ms1_pure_spectrum <- ms1_pure_spectrum[, c('mz_value', 'intensity')]
  ms1_pure_spectrum <- ms1_pure_spectrum[ms1_pure_spectrum['intensity'] != 0, ]
  
  return(list(
    'ms1_pure_spectrum' = ms1_pure_spectrum,
    'comp_ms1' = choosen_comp_ms1  ))
}


#' Select the good MS2 spectrum related to the MS1 peak.
#'
#' @param W_ms2 MS2 pure spectra `matrix`.
#' @param choosen_comp `numeric(1)` good component index.
#'
#' @return Pure MS2 spectrum (normalized) `data.frame`.
#' 
#' @export
#' 
#' @examples
#'  W <- matrix(c(0.2, 5.0, 0.6, 0.4, 0.1, 7.0, 0.3, 0.5, 0.3, 6.0, 0.4, 0.2), nrow = 4, byrow = TRUE)
#'  rownames(W) <- c(1,2,3,4)
#'  comp <- 2
#'  
#'  choose_ms2_pure_spectrum(W_ms2 = W, choosen_comp = comp)
#' 
#' @import magrittr
choose_ms2_pure_spectrum <- function(W_ms2, choosen_comp){
  
  ms2_pure_spectra <- prepare_pure_spectra(W_ms2);
  
  ms2_pure_spectra <-  ms2_pure_spectra %>%
    group_by(comp_nb) 
  # %>%
  #   mutate(intensity = intensity / max(intensity));
  
  ms2_pure_spectrum <- subset(ms2_pure_spectra, comp_nb == paste0('comp',choosen_comp));
  ms2_pure_spectrum <- ms2_pure_spectrum[, c('mz_value','intensity')];
  
  return(ms2_pure_spectrum)
}