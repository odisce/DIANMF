# this function is not completely finished. Improvements are in process!!!!!!!!!!!!!!!!!

#' main function; input: some parameters and the raw data of the mzML file, output: list of features
#'
#' @param rawData.onDiskMSnExp 
#' @param ms1_peaks.df 
#' @param ppm 
#' @param rt_index 
#' @param nmf_parameters.l 
#' @param rank.method 
#' @param initialization.method 
#' @param rt_tol 
#'
#' @return list of identified features
#' @export
dia.nmf.f <- function( # general parameters 
  rawData.onDiskMSnExp = rawData.onDiskMSnExp, ms1_peaks.df = eics_peaks.df,
  # parameters to extract MS1 & MS2 mixed matrices
  ppm = 7, rt_index = TRUE, 
  # NMF parameters
  # to do: add different NMF algorithms 
  nmf_parameters.l = nmf_parameters.l, rank.method = 'concordance', initialization.method = "nndsvd",
  # rt tolerance is used to find MS1 ions that also detect as MS1 peaks
  rt_tol = 0.1){
  
  ms1_peaks.df$is_ion <- FALSE
  peak.idx <- 1
  k <- 1
  features.l <- list()
  
  while( peak.idx <= nrow(ms1_peaks.df) ){
    if( isFALSE(ms1_peaks.df[peak.idx, 'is_ion']) ){
      
      print(peak.idx)
      ms1_mat <- extract_ms_matrix.f(idx.pg = peak.idx, eics_peaks.mat = ms1_peaks.df, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                     ppm = ppm, rt_index = TRUE, mz_range = NULL, iso_win_index = NULL)
      # Delete the rows that have any 1 non-zero value.
      rows_to_keep <- apply(ms1_mat, 1, function(row) sum(row != 0)) > 1
      ms1_mat <- ms1_mat[rows_to_keep, ]
      
      ms1_rank <- 3 # to be modified later on, it should be detected automatically
      ngmcas_res <- nGMCAs(originalMatrix = ms1_mat, rank = ms1_rank, options = nmf_parameters.l,
                           errors_print = nmf_parameters.l$convergence_errors_print, initialization_method = 'nndsvd')
      W_ms1 <- ngmcas_res$S
      H_ms1 <- ngmcas_res$A
      
      # extract pure MS1 spectra and choose the peak corresponding spectrum
      mz_prec <- ms1_peaks.df[peak.idx, 'mz']
      rt_prec <- ms1_peaks.df[peak.idx, 'rt']
      
      mz_ms1_ions <- round(as.numeric(rownames(W_ms1)), 4);
      rownames(W_ms1) <- NULL
      ms1_pure_spectra <- prepare_pure_spectra(W_ms1 , mz_ms1_ions)
      
      # choose the component
      closest_row <- which.min(abs(mz_ms1_ions - mz_prec));
      choosen_comp <- which.max(W_ms1[closest_row, ]);
      ms1_pure_spectrum <- ms1_pure_spectra[ms1_pure_spectra$comp_nb == paste0('comp',choosen_comp), ]
      ms1_pure_spectrum$intensity <- ms1_pure_spectrum$intensity / max(ms1_pure_spectrum$intensity)
      ms1_pure_spectrum <- ms1_pure_spectrum[, c('mz_value', 'intensity')]
      
      # extract MS2 matrices
      min_mz <- min(ms1_pure_spectrum$mz_value)
      max_mz <- max(ms1_pure_spectrum$mz_value)
      info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
      idx.swath <- which(info.swath$lowerMz %between% c(min_mz, max_mz) | info.swath$upperMz %between% c(min_mz, max_mz)) # check this for other values
      
      # so from every isolation window in idx.swath we will extract the ms2 data
      res <- lapply(idx.swath, function(i){
        # i <- 4
        ms2_mat <- extract_ms_matrix.f(idx.pg = peak.idx, eics_peaks.mat = ms1_peaks.df, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                       ppm = ppm, rt_index = TRUE, mz_range = NULL, iso_win_index = i)
        return(ms2_mat)
      })
      
      # Delete the rows that have any 1 non-zero value.
      res2 <- lapply(res, filter_rows)
      
      # apply NMF on all MS1 and MS2 data
      # concatenate MS1 and MS2 matrices
      ms2_matrices <- do.call(rbind, res)
      all_matrices <- rbind(ms1_mat, ms2_matrices)
      
      # I should change the algorithm, add the initialization method!!!!!!!!!!!!!!!!!!
      ms2_rank <- 3  # change this !!!!!!!!!
      rownames(H_ms1) <- NULL
      ngmcas_res_all <- nGMCAs(originalMatrix = all_matrices, rank = ms2_rank, options = nmf_parameters.l,
                               errors_print = nmf_parameters.l$convergence_errors_print, initialization_method = 'subSample', H_sub = H_ms1)
      W <- ngmcas_res_all$S
      H <- ngmcas_res_all$A
      
      # extract pure MS2 spectra and choose the peak corresponding spectrum
      mz_ms2_ions <- round(as.numeric(rownames(ms2_matrices)),4)
      mz_values <- c(mz_ms1_ions, mz_ms2_ions);
      
      closest_row <- which.min(abs(mz_values - mz_prec));
      choosen_comp <- which.max(W[closest_row, ]);   # I think the component is the same as selected in MS1 step
      
      rownames(W) <- NULL
      ms2_pure_spectra <- prepare_pure_spectra(W, mz_values)
      
      ms2_pure_spectra <-  ms2_pure_spectra %>%
        group_by(comp_nb) %>%
        slice(-(1:nrow(ms1_mat))) %>%
        mutate(intensity = intensity / max(intensity)) 
      
      ms2_pure_spectrum <- subset(ms2_pure_spectra, comp_nb == paste0('comp',choosen_comp));
      ms2_pure_spectrum <- ms2_pure_spectrum[, c('mz_value','intensity')]
      
      # so now I have for the peak its pure MS1 and MS2 spectra
      feature_sub.l <- list(
        'peak' = ms1_peaks.df[peak.idx, ],
        'MS1_pure_spectrum' = ms1_pure_spectrum,
        'MS2_pure_spectrum' = ms2_pure_spectrum
      )
      features.l[[k]] <- feature_sub.l
      
      # The other components should be also extracted !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      # update the matrix of MS1 peaks
      ions_peak <- check_ms1_ions(spectrum.df = ms1_pure_spectrum, peaks.mat = ms1_peaks.df, rt_prec = rt_prec, rt_tol = rt_tol)
      # these ions will be deleted from the MS1 peaks matrix, so they are not factorized again, but they may be used in different peaks factorization
      ms1_peaks.df[ions_peak, 'is_ion'] <- TRUE 
      
      k <- k +1
      
    }else{   # peaks which are ions will not be factorized alone
      peak.idx <- peak.idx + 1
      next
    }
    
    peak.idx <- peak.idx + 1
  }
  
  return(features.l)
}