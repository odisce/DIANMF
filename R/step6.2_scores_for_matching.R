# this script contains the 3 classical scores: dot product, inverse dot product and the % of fragments (calculated as described in MS-Dial)

#' presence of fragments score: calculate the percentage of library spectrum fragments found in the experimental spectrum
#'
#' @param measuredSpectra list experimental spectrum
#' @param librarySpectra  list reference spectrum
#' @param bin numeric mz tolerance
#'
#' @return numeric percentage of library fragments in the experimental spectrum
#'
GetPresenceSimilarity <- function(measuredSpectra, librarySpectra, bin){

  if(nrow(librarySpectra) == 0) {return(0)}
  sumM <- sumL <- 0
  minMz <- librarySpectra[1,'mz_value']
  maxMz <- librarySpectra[nrow(librarySpectra), 'mz_value']
  focusedMz <- minMz
  maxLibIntensity <- max(librarySpectra$intensity)
  remaindIndexM <- remaindIndexL <- counter <- libCounter <- 0

  measuredMassList <- list()
  referenceMassList <- list()

  while(focusedMz <= maxMz){
    sumL = 0
    for(i in 1:nrow(librarySpectra)){
      if(librarySpectra[i, 'mz_value'] < focusedMz - bin) {next}
      else if(focusedMz - bin <= librarySpectra[i, 'mz_value'] && librarySpectra[i, 'mz_value'] < focusedMz + bin)
        sumL = sumL + librarySpectra[i, 'intensity']
      else {
        remaindIndexL = i
        break }   }
    if(sumL >= 0.01 * maxLibIntensity) {
      libCounter = libCounter + 1   }

    sumM = 0
    for(i in 1:nrow(measuredSpectra)) {
      if(measuredSpectra[i, 'mz_value'] < focusedMz - bin) { next }
      else if(focusedMz - bin <= measuredSpectra[i, 'mz_value'] && measuredSpectra[i, 'mz_value'] < focusedMz + bin){
        sumM = sumM + measuredSpectra[i, 'intensity']
      } else {
        remaindIndexM = i
        break }  }

    if(sumM > 0 && sumL >= 0.01 * maxLibIntensity) {
      counter = counter +1 }

    if(focusedMz + bin > librarySpectra[nrow(librarySpectra), 'mz_value']) { break }
    focusedMz = librarySpectra[remaindIndexL, 'mz_value']   }

  if(libCounter == 0 ) { return(0) }
  else
  { return(counter / libCounter) }
}

#-------------------------------------------------------------------------------

#' dot product score
#'
#' @param measuredSpectra list experimental spectrum
#' @param librarySpectra  list reference spectrum
#' @param bin numeric mz tolerance
#'
#' @return numeric dot product
#'
GetSimpleDotProductSimilarity <- function(measuredSpectra, librarySpectra, bin){
  scalarM <- scalarR <- covariance <- sumM <-  sumR <- 0

  if(nrow(measuredSpectra) == 0) { return(0) }
  if(nrow(librarySpectra) == 0) { return(0) }

  minMz <- min( measuredSpectra[1, 'mz_value'] ,librarySpectra[1,'mz_value'] )
  maxMz <- max( measuredSpectra[nrow(measuredSpectra), 'mz_value'], librarySpectra[nrow(librarySpectra), 'mz_value'])

  focusedMz <- minMz
  remaindIndexM <- remaindIndexL <- 0

  measuredMassList <- list()
  referenceMassList <- list()

  sumMeasure <- sumReference <- 0
  baseM <- baseR <- -Inf

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
      break }
    }

    if( sumM <= 0 && sumR > 0){
      measuredMassList <- c( measuredMassList, list(c(focusedMz, sumM)) )
      if (sumM > baseM) { baseM = sumM }

      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumR)))
      if(sumR > baseR) { baseR = sumR }
    } else {
      measuredMassList <- c( measuredMassList, list(c(focusedMz, sumM)) )
      if (sumM > baseM) { baseM = sumM }

      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumR)))
      if(sumR > baseR) { baseR = sumR }
    }

    if(focusedMz + bin > max(measuredSpectra[nrow(measuredSpectra) ,'mz_value'], librarySpectra[nrow(librarySpectra), 'mz_value']))
    { break }
    if(focusedMz + bin > librarySpectra[remaindIndexL, 'mz_value'] && focusedMz + bin <= measuredSpectra[remaindIndexM, 'mz_value'])
    { focusedMz = measuredSpectra[remaindIndexM, 'mz_value']
    } else if(focusedMz + bin <= librarySpectra[remaindIndexL, 'mz_value'] && focusedMz + bin > measuredSpectra[remaindIndexM, 'mz_value']){
      focusedMz = librarySpectra[remaindIndexL, 'mz_value']
    } else {
      focusedMz = min(measuredSpectra[remaindIndexM, 'mz_value'], librarySpectra[remaindIndexL, 'mz_value'])}

  }

  if (baseM == 0 || baseR == 0) { return(0) }

  eSpectrumCounter = 0
  lSpectrumCounter = 0

  for(i in 1:length(measuredMassList)) {
    measuredMassList[[i]][2] = measuredMassList[[i]][2] / baseM
    referenceMassList[[i]][2] = referenceMassList[[i]][2] / baseR
    sumMeasure = sumMeasure + measuredMassList[[i]][2]
    sumReference = sumReference + referenceMassList[[i]][2]

    if (measuredMassList[[i]][2] > 0.1) { eSpectrumCounter = eSpectrumCounter + 1 }
    if (referenceMassList[[i]][2] > 0.1){ lSpectrumCounter = lSpectrumCounter + 1}
  }

  peakCountPenalty = 1.0
  if (lSpectrumCounter == 1) {peakCountPenalty = 0.75}
  else if (lSpectrumCounter == 2) {peakCountPenalty = 0.88}
  else if (lSpectrumCounter == 3) {peakCountPenalty = 0.94}
  else if (lSpectrumCounter == 4) {peakCountPenalty = 0.97}

  wM <- wR <- 0
  if (sumMeasure - 0.5 == 0) {wM = 0} else {wM = 1 / (sumMeasure - 0.5)}

  if (sumReference - 0.5 == 0) {wR = 0} else {wR = 1 / (sumReference - 0.5)}

  cutoff = 0.01
  for (i in 1:length(measuredMassList)) {
    if (measuredMassList[[i]][2] < cutoff) { next }

    scalarM = scalarM + measuredMassList[[i]][2] * measuredMassList[[i]][1];
    scalarR = scalarR + referenceMassList[[i]][2] * referenceMassList[[i]][1];
    covariance = covariance + sqrt(measuredMassList[[i]][2] * referenceMassList[[i]][2]) * measuredMassList[[i]][1];
  }

  if(scalarM == 0 || scalarR == 0) { return(0) }
  else {
    return( (covariance^2) / scalarM / scalarR * peakCountPenalty )
  }
}

#-------------------------------------------------------------------------------

#' inverse dot product score
#'
#' @param measuredSpectra list experimental spectrum
#' @param librarySpectra  list reference spectrum
#' @param bin numeric mz tolerance
#'
#' @return numeric inverse dot product
#'
getReverseSearchingSimilarity <- function(measuredSpectra, librarySpectra, bin){

  scalarM <- scalarR <- covariance <- sumM <- sumL <- 0

  if(nrow(librarySpectra) == 0) return(0)

  minMz <- librarySpectra[1,'mz_value']
  maxMz <- librarySpectra[nrow(librarySpectra), 'mz_value']

  focusedMz = minMz; remaindIndexM <- remaindIndexL <- counter <- 0;
  measuredMassList <- list()
  referenceMassList <- list()

  sumMeasure <- sumReference <- 0; baseM <- baseR <- -Inf;

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

    if(sumM <= 0) {
      measuredMassList <- c(measuredMassList, list(c(focusedMz, sumM)))
      if(sumM > baseM) { baseM = sumM }

      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumL)))
      if(sumL > baseR){ baseR = sumL }
    } else {

      measuredMassList <- c(measuredMassList, list(c(focusedMz, sumM)))
      if(sumM > baseM){ baseM <- sumM }

      referenceMassList <- c(referenceMassList, list(c(focusedMz, sumL)))
      if(sumL > baseR) { baseR <- sumL }

      counter = counter + 1
    }

    if(focusedMz + bin > librarySpectra[nrow(librarySpectra), 'mz_value']) { break }
    focusedMz <- librarySpectra[remaindIndexL, 'mz_value']
  }

  if(baseM == 0 || baseR == 0) { return(0) }

  eSpectrumCounter <- lSpectrumCounter <- 0
  for(i in 1: length(measuredMassList)) {
    measuredMassList[[i]][2] = measuredMassList[[i]][2] / baseM
    referenceMassList[[i]][2] = referenceMassList[[i]][2] / baseR
    sumMeasure = sumMeasure + measuredMassList[[i]][2]
    sumReference = sumReference + referenceMassList[[i]][2]

    if(measuredMassList[[i]][2] > 0.1) {eSpectrumCounter = eSpectrumCounter + 1}
    if(referenceMassList[[i]][2] > 0.1) {lSpectrumCounter = lSpectrumCounter + 1}
  }

  peakCountPenalty  <- 1.0
  if(lSpectrumCounter == 1){
    peakCountPenalty = 0.75
  } else if(lSpectrumCounter == 2){
    peakCountPenalty = 0.88
  } else if(lSpectrumCounter == 3){
    peakCountPenalty = 0.94
  } else if(lSpectrumCounter == 4){
    peakCountPenalty = 0.97 }

  wM <- wR <- 0

  if(sumMeasure - 0.5 == 0){
    wM = 0} else{
      wM = 1 / (sumMeasure - 0.5) }

  if(sumReference - 0.5 == 0){
    wR = 0 } else {
      wR = 1 / (sumReference - 0.5)}

  cutoff = 0.01

  for(i in 1:length(measuredMassList)) {
    if(referenceMassList[[i]][2] < cutoff){
      next }

    scalarM = scalarM + measuredMassList[[i]][2] * measuredMassList[[i]][1]
    scalarR = scalarR + referenceMassList[[i]][2] * referenceMassList[[i]][1]
    covariance = covariance + sqrt(measuredMassList[[i]][2] * referenceMassList[[i]][2]) * measuredMassList[[i]][1]

    # scalarM = scalarM + measuredMassList[[i]][2]
    # scalarR = scalarR + referenceMassList[[i]][2]
    # covariance = covariance + sqrt(measuredMassList[[i]][2] * referenceMassList[[i]][2])
  }

  if(scalarM == 0 || scalarR == 0) { return(0) }
  else { return( (covariance^2) / scalarM / scalarR * peakCountPenalty) }
}

#-------------------------------------------------------------------------------

# match_pure_scores2 <- function(polarity, mz_precursor, data_base, measured_spectra){
#
#   mz_tol <- mz_precursor * 5e-6
#   idx_ref_spect <- data_base$index[Polarity == polarity & dplyr::between(mz, mz_precursor-mz_tol, mz_precursor+mz_tol), id]
#
#   if( length(idx_ref_spect) > 0 ){
#     res <- lapply(seq_along(idx_ref_spect), function(i){
#
#       print(paste('spec index', i, '-->', idx_ref_spect[i]))
#
#       ref_spect <- data_base$spectra[[ idx_ref_spect[i] ]]@spectrum
#       # # normalize the spectrums
#       # measured_spectra$intensity <- measured_spectra$intensity / max(measured_spectra$intensity)
#
#       ref_spect <- as.data.frame(ref_spect)
#       colnames(ref_spect) <- c('mz_value', 'intensity')
#       ref_spect$intensity <- ref_spect$intensity / max(ref_spect$intensity)
#
#       score1 <- GetSimpleDotProductSimilarity(measured_spectra, ref_spect,  bin = mz_tol)
#       score2 <- GetPresenceSimilarity(measured_spectra, ref_spect, bin = mz_tol)
#       score3 <- getReverseSearchingSimilarity(measured_spectra, ref_spect, bin = mz_tol)
#       score = (score1 + score2 + score3) / 3
#
#       data.frame(
#         "ref spectrum index" = idx_ref_spect[i],
#         "mz" = data_base$index[idx_ref_spect[i], 'mz'],
#         # "Polarity" = polarity,
#         "NCE" =  data_base$index[idx_ref_spect[i], 'NCE'],
#         "Composition" = data_base$index[idx_ref_spect[i], 'composition'],
#         "name" = data_base$index[idx_ref_spect[i], 'name'],
#         "Score1" = score1,
#         "Score2" = score2,
#         "Score3" = score3,
#         "Total Score" = score)
#     })
#     res1 <- do.call(rbind, res)
#   } else { return(NULL) }
#   return(res1)
# }

# --------------- trying to make the running time faster -----------------------

#' match peaks with some reference spectra from Newdb.rds
#'
#' @param polarity string 'POS' or 'NEG'
#' @param mz_precursor numeric peak/MS1 precursor mz value
#' @param data_base reference/library database
#' @param measured_spectra list experimental pure spectrum
#' @param mz_tol numeric
#'
#' @importFrom dplyr between
#'
#' @return list matching scores with some information of the matched library(reference) spectrum
#'
match_pure_scores2 <- function(polarity, mz_precursor, data_base, measured_spectra, mz_tol = NULL) {

  if(is.null(mz_tol)){ # automatic mz_tol depending on the precursor mz
    mz_tol <- mz_precursor * 5e-6
  } # else the mz_tol given by the user
  
  idx_ref_spect <- data_base$index[Polarity == polarity & between(mz, mz_precursor-mz_tol, mz_precursor+mz_tol), id]
  # print(idx_ref_spect)
  
  if (length(idx_ref_spect) > 0) {
    # # Normalize measured spectra intensities
    # measured_spectra$intensity <- measured_spectra$intensity / max(measured_spectra$intensity)

    results <- mclapply(idx_ref_spect, function(i) {
      cat("Processing ref_spect index:", i, "\n")

      ref_spect <- data_base$spectra[[i]]@spectrum
      ref_spect <- as.data.frame(ref_spect)
      colnames(ref_spect) <- c('mz_value', 'intensity')
      # Normalize reference spectra intensities
      ref_spect$intensity <- ref_spect$intensity / max(ref_spect$intensity)

      score1 <- GetSimpleDotProductSimilarity(measured_spectra, ref_spect, bin = mz_tol)
      score2 <- GetPresenceSimilarity(measured_spectra, ref_spect, bin = mz_tol)
      score3 <- getReverseSearchingSimilarity(measured_spectra, ref_spect, bin = mz_tol)
      # score <- ( score1 + score2 + score3 ) /3
      score <- ( score2 + score3 ) /2

      data.table(
        "ref spectrum index" = i,
        "mz" = data_base$index[i, 'mz'],
        "NCE" = data_base$index[i, 'NCE'],
        "Composition" = data_base$index[i, 'composition'],
        "name" = data_base$index[i, 'name'],
        "Dot prod" = round(score1,2),
        "Rev prod" = round(score3,2),
        "Presence" = round(score2,2),
        "total score" = round(score,2) )
    })
    res1 <- rbindlist(results)
  } else {
    return(NULL)
  }
  return(res1)
}

# match with the whole data base
# match_with_whole_data_base <- function(polarity, mz_precursor, data_base, measured_spectra) {
#
#   mz_tol <- mz_precursor * 5e-6
#   # mz_tol <- 0.01
#   idx_ref_spect <- data_base$index[Polarity == polarity, id]
#
#   if (length(idx_ref_spect) > 0) {
#     # # Normalize measured spectra intensities
#     # measured_spectra$intensity <- measured_spectra$intensity / max(measured_spectra$intensity)
#
#     results <- mclapply(idx_ref_spect, function(i) {
#       cat("Processing ref_spect index:", i, "\n")
#
#       ref_spect <- data_base$spectra[[i]]@spectrum
#       ref_spect <- as.data.frame(ref_spect)
#       colnames(ref_spect) <- c('mz_value', 'intensity')
#       # Normalize reference spectra intensities
#       ref_spect$intensity <- ref_spect$intensity / max(ref_spect$intensity)
#
#       score1 <- GetSimpleDotProductSimilarity(measured_spectra, ref_spect, bin = mz_tol)
#       score2 <- GetPresenceSimilarity(measured_spectra, ref_spect, bin = mz_tol)
#       score3 <- getReverseSearchingSimilarity(measured_spectra, ref_spect, bin = mz_tol)
#       score <- (score1 + score2 + score3) /3
#
#       data.table(
#         "ref spectrum index" = i,
#         "mz" = data_base$index[i, 'mz'],
#         "NCE" = data_base$index[i, 'NCE'],
#         "Composition" = data_base$index[i, 'composition'],
#         "name" = data_base$index[i, 'name'],
#         "Dot prod" = round(score1,3),
#         "Rev prod" = round(score3,3),
#         "Presence" = round(score2,3),
#         "total score" = round(score,3) )
#     })
#     res1 <- rbindlist(results)
#   } else {
#     return(NULL)
#   }
#   return(res1)
# }

#-------------------------------------------------------------------------------
