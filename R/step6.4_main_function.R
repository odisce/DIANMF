# this function is not completely finished. Improvements are in process!!!!!!!!!!!!!!!!!

#' main function; input: some parameters and the raw data of the mzML file, output: list of features
#'
#' @param mzML_path character mzML file path
#' @param MS_level character MS1 or MS2
#' @param ms1_peaks matrix or data.frame contains all MS1 peaks information
#' @param ppm numeric(1)
#' @param peaks_nb numeric
#' @param rt_index Boolean
#' @param maximumIteration numeric
#' @param maxFBIteration numeric 
#' @param toleranceFB numeric
#' @param MS1_init_method character
#' @param MS2_init_method character
#' @param errors_print Boolean
#' @param rt_tol character
#' @param plot_key Boolean
#' @param plot_path character
#'
#' @return list of identified features
#' @export
#' @importFrom purrr map
#' @import magrittr
#' @import dplyr
#' @importFrom data.table %between%
#' @importFrom MSnbase readMSData
dia.nmf.f <- function(
    mzML_path = NULL,
    MS_level = "MS2", ppm = 7,
    # parameters to detect peaks by xcms, or input the peak matrix or data.frame
    ms1_peaks = NULL,
    # parameter to identify all MS1 peaks or specific nb
    peaks_nb = NULL,
    # parameters to extract MS1 & MS2 mixed matrices
    rt_index = TRUE, 
    # NMF parameters
    maximumIteration = 10, maxFBIteration = 5, toleranceFB = 1e-5,
    MS1_init_method = 'nndsvd', MS2_init_method = 'subSample', errors_print = FALSE,
    # rt tolerance is used to find MS1 ions that also detect as MS1 peaks
    rt_tol = 2,
    # to plot all mixed and pure MS1 & MS2 data
    plot_key = FALSE, plot_path = NULL 
    ){
  
  rawData.onDiskMSnExp <- MSnbase::readMSData(mzML_path, mode = "onDisk");
  ms1_peaks.df <- prepare_ms1_peaks(ms1_peaks = ms1_peaks)

  peak.idx <- 1
  k <- 1
  features.l <- list()
  
  if( is.null(peaks_nb) ){
    peaks_to_be_identified <-  nrow(ms1_peaks.df)
  } else {
    peaks_to_be_identified <- peaks_nb
  };
  
  while( peak.idx <= peaks_to_be_identified ){
    if( isFALSE(ms1_peaks.df[peak.idx, 'is_ion']) ){
      
      # print(peak.idx)
      
      ms1_mat <- extract_ms_matrix.f(idx.pg = peak.idx, eics_peaks = ms1_peaks.df, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                     ppm = ppm, rt_index = TRUE, mz_range = NULL, iso_win_index = NULL)
      if( nrow(ms1_mat) == 0 ){
        print('No MS1 data.')
        peak.idx <- peak.idx + 1
        next
      };
      
      # NMF on MS1 data
      ms1_rank <- 3 # to be modified later on, it should be detected automatically
      ngmcas_res <- nGMCAs(originalMatrix = ms1_mat, rank = ms1_rank,
                           maximumIteration = maximumIteration, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB,
                           initialization_method = MS1_init_method,
                           errors_print = errors_print);
      W_ms1 <- ngmcas_res$S
      H_ms1 <- ngmcas_res$A
      W_ms1 <- apply(W_ms1, 2, function(x) x / max(x));  # normalize every column (comp spectrum) by is max value
      
      # prepare pure MS1 spectra and choose the peak corresponding spectrum (the good spectrum)
      mz_prec <- ms1_peaks.df[peak.idx, 'mz'];
      rt_prec <- ms1_peaks.df[peak.idx, 'rt'];
      
      # choose the good MS1 spectrum
      ms1_pure_data <- choose_ms1_pure_spectrum(W_ms1, mz_prec);
      ms1_pure_spectrum <- ms1_pure_data$ms1_pure_spectrum;
      choosen_comp_ms1 <- ms1_pure_data$choosen_comp_ms1;
      
      # Test the ions in the chosen ms1 pure spectra, which will also be considered as peaks or not.-----------------------------------------
      ions_maybe_peaks <- which(W_ms1[,choosen_comp_ms1] >= 0.6 * rowSums(W_ms1) );
      ions_maybe_peaks <- as.numeric(names(ions_maybe_peaks));
      ions_are_peaks <- check_ms1_ions(ions_maybe_peaks.v = ions_maybe_peaks, peaks.df = ms1_peaks.df, rt_prec = rt_prec, rt_tol = rt_tol);
      # these ions are not factorized again, but they may be used in different peaks factorization
      ms1_peaks.df[ions_are_peaks, 'is_ion'] <- TRUE;
      # --------------------------------------------------------------------------------------------------------- the peaks df is updated
      
      if ( MS_level == "MS1" ){
        
        feature_sub.l <- list( 
          'peak' = ms1_peaks.df[peak.idx, ],
          'MS1_pure_spectrum' = ms1_pure_spectrum,
          'MS1_mixed_mat' = ms1_mat,
          'W_ms1' = W_ms1,
          'H_ms1' = H_ms1,
          'comp_ms1' = choosen_comp_ms1,
          'ms1_ions_are_peaks' = ions_are_peaks  );
        
      } else { # the user want to process the MS2 data also
        
        info.swath <- isolationWindows.range(rawData.onDiskMSnExp);
        
        res_ms2 <- extract_ms2_matrices(peak.idx = peak.idx, ms1_peaks.df = ms1_peaks.df,
                                        ppm = ppm, ms1_pure_spectrum = ms1_pure_spectrum,
                                        rawData.onDiskMSnExp = rawData.onDiskMSnExp, info.swath = info.swath )
        ms2_matrices <- do.call(rbind, res_ms2);
        if( nrow(ms2_matrices) == 0 ){
          print('No MS2 data.')
          peak.idx <- peak.idx + 1
          next
        };

        # NMF on MS2 data
        ms2_rank <- 3  # change this !!!!!!!!!
        rownames(H_ms1) <- NULL; 
        ngmcas_res_all <- nGMCAs(originalMatrix = ms2_matrices, rank = ms2_rank,
                                 maximumIteration = maximumIteration, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB,
                                 initialization_method = MS2_init_method, H_sub = H_ms1,
                                 errors_print = errors_print);
        W_ms2 <- ngmcas_res_all$S;
        H_ms2 <- ngmcas_res_all$A;
        W_ms2 <- apply(W_ms2, 2, function(x) x / max(x));  # normalize every column (comp spectrum) by is max value
        
        # choose the good MS2 spectrum
        choosen_comp_ms2 <- elutions_corelation(chromo_main = H_ms1[choosen_comp_ms1, ], chromos = H_ms2);
        ms2_pure_spectrum <- choose_ms2_pure_spectrum(W_ms2, choosen_comp_ms2);
        
        # Filter the ms2 pure spectra to get only the fragments related to the precursor from the SWATH window where it was fragmented
        ms2_pure_spectrum_specific <- filter_ms2_spectrum(ms2_pure_spectrum = ms2_pure_spectrum, ms2_matrices = res_ms2, mz_prec, info.swath = info.swath );
        if( !is.null(ms2_pure_spectrum_specific) ){
          ms2_pure_spectrum_specific <- ms2_pure_spectrum_specific[ms2_pure_spectrum_specific['intensity'] != 0, ];
        };
        
        # now I can delete the zero intensity fragments
        ms2_pure_spectrum <- ms2_pure_spectrum[ms2_pure_spectrum['intensity'] != 0, ];
        
        if( isTRUE(plot_key) & !is.null(plot_path) & isTRUE(rt_index) ){
          plot_all_info(d.plot = plot_path,
                        peak.idx = peak.idx,
                        ms_mixed1 = ms1_mat, ms_mixed2 = ms2_matrices,
                        ms_pure_H1 = H_ms1, ms_pure_H2 = H_ms2,
                        ms_pure_W1 = W_ms1, ms_pure_W2 = W_ms2,
                        rt_prec = rt_prec, mz_prec = mz_prec,
                        choosen_comp_ms1 = choosen_comp_ms1, choosen_comp_ms2 = choosen_comp_ms2);
        }
        
        # save the peak info, this part it will be deleted at the end
        feature_sub.l <- list(
          'peak' = ms1_peaks.df[peak.idx, ],
          'MS1_pure_spectrum' = ms1_pure_spectrum,
          'MS2_pure_spectrum' = ms2_pure_spectrum,
          'MS2_pure_spectrum_specific' = ms2_pure_spectrum_specific,
          'MS1_mixed_mat' = ms1_mat,
          'MS2_mixed_mat' = res_ms2,
          'W_ms1' = W_ms1,
          'H_ms1' = H_ms1,
          'W_ms2' = W_ms2,
          'H_ms2' = H_ms2,
          'comp_ms1' = choosen_comp_ms1,
          'comp_ms2' = choosen_comp_ms2,
          'ions_are_peaks' = ions_are_peaks  );
        
      }
      
      # the following line will be deleted at the end !!!!!!!!!!!!!!!!!!!!!!!!!!
      # saveRDS(feature_sub.l, paste0('C:/Users/DK273056/Documents/DIA_NMF_R_package_outputs_new/peak_', peak.idx, '.rds') );
      
      features.l[[k]] <- feature_sub.l;
      k <- k + 1
      } else {   # peaks which are ions will not be factorized alone
          peak.idx <- peak.idx + 1
          next
        }
    
      peak.idx <- peak.idx + 1
  }
  
  return(features.l)
}
