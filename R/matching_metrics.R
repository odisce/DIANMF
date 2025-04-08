# classical scores

#' Calculate the Dot product
#'
#' @param measuredSpectra experimental spectrum.
#' @param librarySpectra reference spectrum.
#' @param bin `numeric` tolerance to merge.
#'
#' @returns `numeric` dot product score.
GetSimpleDotProductSimilarity <- function(measuredSpectra, librarySpectra, bin){
  scalarM <- scalarR <- covariance <- sumM <- sumR <- 0
  
  if(nrow(measuredSpectra) == 0) { return(0) }
  if(nrow(librarySpectra) == 0) { return(0) }
  
  measuredSpectra <- measuredSpectra[order(measuredSpectra$mz_value), ]
  librarySpectra <- librarySpectra[order(librarySpectra$mz_value), ]
  
  minMz <- min( measuredSpectra[1, 'mz_value'] ,librarySpectra[1,'mz_value'] )
  maxMz <- max( measuredSpectra[nrow(measuredSpectra), 'mz_value'], librarySpectra[nrow(librarySpectra), 'mz_value'])
  
  focusedMz <- minMz
  remaindIndexM <- remaindIndexL <- 0
  
  measuredMassList <- list()
  referenceMassList <- list()
  
  while(focusedMz <= maxMz) {
    sumM = 0
    for(i in  1:nrow(measuredSpectra)) {
      if(measuredSpectra[i,'mz_value'] < focusedMz - bin) { next }
      else if(focusedMz - bin <= measuredSpectra[i,'mz_value'] && measuredSpectra[i,'mz_value'] < focusedMz + bin){
        sumM = sumM + measuredSpectra[i,'intensity']
      } else {
        remaindIndexM = i
        break }
    }
    
    sumR = 0
    for(i in 1:nrow(librarySpectra)){
      if (librarySpectra[i,'mz_value'] < focusedMz - bin) { next }
      else if(focusedMz - bin <= librarySpectra[i,'mz_value'] && librarySpectra[i,'mz_value'] < focusedMz + bin){
        sumR = sumR + librarySpectra[i,'intensity']
      } else { remaindIndexL = i
      break }  }
    
    # get the max intensity of the peak for normalization
    if( sumM <= 0 && sumR > 0){
      measuredMassList <- c( measuredMassList, list(c(focusedMz, sumM)) )
      
      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumR)))
    } else {
      measuredMassList <- c( measuredMassList, list(c(focusedMz, sumM)) )
      
      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumR)))
    }
    
    if(focusedMz + bin > max(measuredSpectra[nrow(measuredSpectra) ,'mz_value'], librarySpectra[nrow(librarySpectra), 'mz_value']))
    { break }
    if(focusedMz + bin > librarySpectra[remaindIndexL, 'mz_value'] && focusedMz + bin <= measuredSpectra[remaindIndexM, 'mz_value'])
    { focusedMz = measuredSpectra[remaindIndexM, 'mz_value']
    } else if(focusedMz + bin <= librarySpectra[remaindIndexL, 'mz_value'] && focusedMz + bin > measuredSpectra[remaindIndexM, 'mz_value']){
      focusedMz = librarySpectra[remaindIndexL, 'mz_value']
    } else {
      focusedMz = min(measuredSpectra[remaindIndexM, 'mz_value'], librarySpectra[remaindIndexL, 'mz_value']) }
  }
  
  for (i in 1:length(measuredMassList)) {
    scalarM = scalarM + measuredMassList[[i]][2] * measuredMassList[[i]][1];
    scalarR = scalarR + referenceMassList[[i]][2] * referenceMassList[[i]][1];
    covariance = covariance + sqrt(measuredMassList[[i]][2] * referenceMassList[[i]][2]) * measuredMassList[[i]][1];
  }
  
  if(scalarM == 0 || scalarR == 0) { 
    return(0)
  } else { 
    return( (covariance^2) / scalarM / scalarR )
  }
}


#' Calculate the matched fragments ratio
#'
#' @inheritParams GetSimpleDotProductSimilarity
#'
#' @returns `numeric` fragments ratio score.
GetPresenceSimilarity <- function(measuredSpectra, librarySpectra, bin){
  
  if(nrow(librarySpectra) == 0) {return(0)}
  
  # the spectra should be ordered in increasing order of intensity
  measuredSpectra <- measuredSpectra[order(measuredSpectra$mz_value), ]
  librarySpectra <- librarySpectra[order(librarySpectra$mz_value), ]
  
  sumM <- sumL <- 0
  minMz <- librarySpectra[1,'mz_value']
  maxMz <- librarySpectra[nrow(librarySpectra), 'mz_value']
  focusedMz <- minMz
  # maxLibIntensity <- max(librarySpectra$intensity)
  remaindIndexM <- remaindIndexL <- counter <- libCounter <- 0
  
  while(focusedMz <= maxMz){
    sumL = 0
    for(i in 1:nrow(librarySpectra)){
      if(librarySpectra[i, 'mz_value'] < focusedMz - bin) {next}
      else if(focusedMz - bin <= librarySpectra[i, 'mz_value'] && librarySpectra[i, 'mz_value'] < focusedMz + bin)
        sumL = sumL + librarySpectra[i, 'intensity']
      else {
        remaindIndexL = i
        break }   }
    
    if(sumL > 0) {
      libCounter = libCounter + 1 }
    
    sumM = 0
    for(i in 1:nrow(measuredSpectra)) {
      if(measuredSpectra[i, 'mz_value'] < focusedMz - bin) { next }
      else if(focusedMz - bin <= measuredSpectra[i, 'mz_value'] && measuredSpectra[i, 'mz_value'] < focusedMz + bin){
        sumM = sumM + measuredSpectra[i, 'intensity']
      } else {
        remaindIndexM = i
        break }  }
    
    if(sumM > 0 && sumL > 0) {
      counter = counter + 1 }
    
    if(focusedMz + bin > librarySpectra[nrow(librarySpectra), 'mz_value']) { break }
    focusedMz = librarySpectra[remaindIndexL, 'mz_value']   }
  
  if(libCounter == 0 ) { 
    return(0)
  } else { 
    return(counter / libCounter)
  }
}


#' Calculate reverse dot-product
#'
#' @inheritParams GetSimpleDotProductSimilarity
#'
#' @returns `numeric` reverse dot product score.
getReverseSearchingSimilarity <- function(measuredSpectra, librarySpectra, bin){
  
  scalarM <- scalarR <- covariance <- sumM <- sumL <- 0
  
  if(nrow(librarySpectra) == 0) return(0)
  
  measuredSpectra <- measuredSpectra[order(measuredSpectra$mz_value), ]
  librarySpectra <- librarySpectra[order(librarySpectra$mz_value), ]
  
  minMz <- librarySpectra[1,'mz_value']
  maxMz <- librarySpectra[nrow(librarySpectra), 'mz_value']  # just against fragments in the library spectrum (not the common fragments)
  
  focusedMz = minMz
  remaindIndexM <- remaindIndexL <- 0
  measuredMassList <- list()
  referenceMassList <- list()
  sumMeasure <- sumReference <- 0
  
  while(focusedMz <= maxMz){
    
    sumM = 0
    for(i in 1:nrow(measuredSpectra)){
      if(measuredSpectra[i, 'mz_value'] < (focusedMz - bin)){
        next
      } else if( (focusedMz - bin) <= measuredSpectra[i, 'mz_value'] && measuredSpectra[i, 'mz_value'] < (focusedMz + bin) ){
        sumM = sumM + measuredSpectra[i, 'intensity']
      } else {
        remaindIndexM = i
        break}
    }
    
    sumL = 0
    for(i in 1:nrow(librarySpectra)){
      if(librarySpectra[i, 'mz_value'] < (focusedMz - bin)){
        next
      } else if( (focusedMz - bin) <= librarySpectra[i, 'mz_value'] && librarySpectra[i, 'mz_value'] < (focusedMz + bin)){
        sumL = sumL + librarySpectra[i, 'intensity']
      } else {
        remaindIndexL = i
        break }
    }
    
    if( sumM > 0 ){
      measuredMassList <- c(measuredMassList, list(c(focusedMz, sumM)))
      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumL))) 
    }
    
    if(focusedMz + bin > librarySpectra[nrow(librarySpectra), 'mz_value']) { break }
    focusedMz <- librarySpectra[remaindIndexL, 'mz_value']
  }
  
  if( length(measuredMassList) == 0 ){
    return(0)  }
  
  for(i in 1:length(measuredMassList)) {
    scalarM = scalarM + measuredMassList[[i]][2] * measuredMassList[[i]][1]
    scalarR = scalarR + referenceMassList[[i]][2] * referenceMassList[[i]][1]
    covariance = covariance + sqrt(measuredMassList[[i]][2] * referenceMassList[[i]][2]) * measuredMassList[[i]][1]
  }
  
  if(scalarM == 0 || scalarR == 0) { 
    return(0)
  } else {
    return( (covariance^2) / (scalarM * scalarR) )
  }
}


ms2_matching <- function(polarity, mz_prec, data_base, measured_spectra, mz_tol = NULL, idx = NA) {

  if ( is.null(mz_tol) ){
    mz_tol <- mz_prec * 5e-6
  }

  if( is.na(idx) ){
    idx_ref_spect <- data_base$index[Polarity == polarity & dplyr::between(mz, mz_prec - 0.005, mz_prec + 0.005), id]
  } else {
    idx_ref_spect <- idx
  }

  # change the class of measured_spectra from Spectra to df
  
  if (length(idx_ref_spect) > 0) {
    # Normalize measured spectra intensities
    measured_spectra$intensity <- measured_spectra$intensity / max(measured_spectra$intensity)

    results <- parallel::mclapply(idx_ref_spect, function(i) {
      cat("Processing ref_spect index:", i, "\n")

      ref_spect <- data_base$spectra[[i]]@spectrum
      ref_spect <- as.data.frame(ref_spect)
      colnames(ref_spect) <- c('mz_value', 'intensity')
      ref_spect <- ref_spect[ ref_spect$intensity >= 0.01 *max(ref_spect$intensity), ]

      # Normalize reference spectra intensities
      ref_spect$intensity <- ref_spect$intensity / max(ref_spect$intensity)
      if( nrow(ref_spect) == 1 ){
        return(NULL)
      }
      score1 <- GetSimpleDotProductSimilarity(measured_spectra, ref_spect, bin = mz_tol)
      score2 <- getReverseSearchingSimilarity(measured_spectra, ref_spect, bin = mz_tol)
      score3 <- GetPresenceSimilarity(measured_spectra, ref_spect, bin = mz_tol)
      score <- ( score2 + score3 ) /2

      data.table(
        "ref_spectrum_index" = i,
        "mz" = as.numeric(data_base$index[i, 'mz']),
        "NCE" = as.character(unname(data_base$index[i, 'NCE'])),
        "Composition" = as.character(unname(data_base$index[i, 'composition'])),
        "name" = as.character(unname(data_base$index[i, 'name'])),
        "dot_prod" = round(score1,3),
        "rev_prod" = round(score2,3),
        "presence" = round(score3,3),
        "total_score" = round(score,3) )
    })
    res1 <- rbindlist(results)
  } else {
    return(NULL)
  }
  return(res1)
}


prepare_spectrum <- function(spect_df){
  spect_df <- spect_df[, c('mz', 'value')]
  colnames(spect_df) <- c("mz_value", "intensity")
  spect_df <- spect_df[spect_df$intensity > 0, ]
  spect_df$intensity <- spect_df$intensity / max(spect_df$intensity)
  spect_df <- as.data.frame(spect_df)
  return(spect_df)
}


#'  Plot measured (mixed and pure) vs. library spectra of a feature
#'
#' @param measuredSpectra_pure `data.frame`
#' @param measuredSpectra_mixed `data.frame`
#' @inheritParams GetSimpleDotProductSimilarity
#' @param scores_pure `data.frame`
#' @param scores_mixed `data.frame`
#' @param mslevel "MS1" or "MS2
#'
#' @returns ggplot2 plot
#' @export
#' 
#' @import ggplot2
#' @importFrom ggpubr ggarrange
plot_spectra_vs <- function(measuredSpectra_pure, measuredSpectra_mixed, librarySpectra, scores_pure, scores_mixed, mslevel) {
  
  if (mslevel == "MS2") {
    # MS2: color by IsoWin, include name in caption
    p_mixed <- ggplot() + 
      geom_linerange(data = measuredSpectra_mixed, aes(x = mz_value, ymin = 0, ymax = intensity, color = as.factor(IsoWin))) +
      geom_linerange(data = librarySpectra, aes(x = mz_value, ymin = -intensity, ymax = 0), color = 'red') +
      theme_bw(base_size = 14) +
      guides(color = guide_legend(title = "IsoWin")) +
      labs(caption = paste("MS2 Mix vs. lib [", scores_pure$name, "]:",
                           " DP:", scores_mixed$dp,
                           " IDP:", scores_mixed$rdp,
                           " FPP:", scores_mixed$fp))
    
    p_pure <- ggplot() + 
      geom_linerange(data = measuredSpectra_pure, aes(x = mz_value, ymin = 0, ymax = intensity, color = as.factor(IsoWin))) +
      geom_linerange(data = librarySpectra, aes(x = mz_value, ymin = -intensity, ymax = 0), color = 'red') +
      theme_bw(base_size = 14) +
      guides(color = guide_legend(title = "IsoWin")) +
      labs(caption = paste("MS2 Pure vs. lib [", scores_pure$name, "]:",
                           " DP:", scores_pure$dp,
                           " IDP:", scores_pure$rdp,
                           " FPP:", scores_pure$fp))
    
  } else {
    # MS1 or other: no coloring by IsoWin
    p_mixed <- ggplot() + 
      geom_linerange(data = measuredSpectra_mixed, aes(x = mz_value, ymin = 0, ymax = intensity), color = 'blue') +
      geom_linerange(data = librarySpectra, aes(x = mz_value, ymin = -intensity, ymax = 0), color = 'red') +
      theme_bw(base_size = 14) +
      labs(caption = paste(mslevel, "Mix vs. lib:",
                           " DP:", scores_mixed$dp,
                           " IDP:", scores_mixed$rdp,
                           " FPP:", scores_mixed$fp)) +
      guides(color = "none")
    
    p_pure <- ggplot() + 
      geom_linerange(data = measuredSpectra_pure, aes(x = mz_value, ymin = 0, ymax = intensity), color = 'blue') +
      geom_linerange(data = librarySpectra, aes(x = mz_value, ymin = -intensity, ymax = 0), color = 'red') +
      theme_bw(base_size = 14) +
      labs(caption = paste(mslevel, "Pure vs. lib:",
                           " DP:", scores_pure$dp,
                           " IDP:", scores_pure$rdp,
                           " FPP:", scores_pure$fp)) +
      guides(color = "none")
  }
  
  ggpubr::ggarrange(p_mixed, p_pure, nrow = 2)
}


