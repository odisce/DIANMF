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
prepare_mixed_data <- function(ms_mixed, mz_values, rts ){
  
  rownames(ms_mixed) <- make.unique(as.character(mz_values))
  colnames(ms_mixed) <- rts
  ms_mixed <- reshape2::melt(ms_mixed)
  colnames(ms_mixed) <- c('mz_value', 'rt', 'intensity')
  
  return(ms_mixed)
}

#' Retrieve pure sources: elution profiles and there spectra
#'
#' @param W `matrix` basis.
#' @param H `matrix` corresponding coefficient. 
#' @param ms_type `character` c('max', 'mean', 'sum'). Indicate how to compute the spectrum, i.e., how to calculate the intensity of ions/fragments of the same mz.
#'
#' @return `list` of pure sources. Each list contains pure elution profiles and the pure spectrum. 
#' 
#' @export 
#' 
#' @import magrittr
pure_sources.f <- function(W, H, ms_type = c('max', 'mean', 'sum')){
  res <- lapply(1:ncol(W), function(s){
    # s <- 13
    spect <- as.matrix(W[, s], ncol = 1)
    elut <- as.matrix(H[s, ], nrow = 1) 
    source <- spect %*% t(elut)
    colnames(source) <- as.numeric(colnames(H))
    rownames(source) <- as.numeric(rownames(W))
    source <- source[rowSums(source != 0) > 0, ]
    
    if( !is.null(nrow(source)) & nrow(source) > 2 ){
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
    } else {
      return(NULL)
    }
    return(res)
  })
}
