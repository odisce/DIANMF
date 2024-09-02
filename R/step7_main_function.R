#' main function; input: some parameters and mzML file path, output: list of features.
#' 
#' @param mzML_path mzML file path.
#' @inheritParams extract_mixedMat
#' @param peaks_by_xcms `Logical` `TRUE` if the user wants to detect MS1 levels using XCMS, `FALSE` if the user will provide the peaks matrix or data frame to be identified.
#' @inheritParams prepare_ms1_peaks
#' @param d.out `character` file path to save the MS1 peaks matrix if detected by XCMS.
#' @inheritParams nGMCAs
#' @param MS1_init_method 'character' MS1 initialization method.
#' @param MS2_init_method 'character' MS2 initialization method.
#' @inheritParams check_ms1_ions
#' @inheritParams extract_ms_matrix.f
#' @param ... Additional parameters passed to \code{\link[xcms]{CentWaveParam}}.
#'
#' @return `list` of identified features.
#' @export
#' @importFrom purrr map
#' @import magrittr
#' @import dplyr 
#' @importFrom data.table %between%
#' @importFrom MSnbase readMSData
dia_nmf.f <- function(
    mzML_path = NULL,
    ms_level = "MS2",
    # parameters to detect peaks by xcms, or input the peak matrix or data.frame
    peaks_by_xcms = TRUE, ms1_peaks = NULL, d.out = NULL,
    # parameters to extract MS1 & MS2 mixed matrices
    ppm.n = 7,
    # NMF parameters
    maximumIteration = 10, maxFBIteration = 10, toleranceFB = 1e-5,
    MS1_init_method = 'nndsvd', MS2_init_method = 'subSample', errors_print = FALSE,
    # rt tolerance is used to find MS1 ions that also detect as MS1 peaks
    rt_tol = 2,
    # additional parameters from xcms::CentWaveParam to detect peaks
    ... ){
  
  rawData.onDiskMSnExp <- MSnbase::readMSData(mzML_path, mode = "onDisk");
  
  if( ms_level == "MS1" & !(1 %in% unique(MSnbase::msLevel(rawData.onDiskMSnExp))) ){
    print("Error, No MS1 data in this file.")
    return(NULL)
  };
  
  if( ms_level == "MS2" & !(2 %in% unique(MSnbase::msLevel(rawData.onDiskMSnExp))) ){
    print("Error, No MS2 data in this file.")
    return(NULL)
  };
  
  if( peaks_by_xcms == TRUE & is.null(ms1_peaks) ){
    ms1_peaks <- detect_peaks_by_xcms(rawData.onDiskMSnExp, ...)
    
    if (!is.null(d.out)) {     # save the peak matrix
      if (!dir.exists(d.out)) {
        dir.create(d.out)
      }
      saveRDS(ms1_peaks, file = paste0(d.out, '/ms1_peaks.rds'))
    }
    
    if( nrow(ms1_peaks) == 0 ){
      print("No detected MS1 peaks!")
      return(NULL)
    }
  };

  ms1_peaks.df <- prepare_ms1_peaks(ms1_peaks = ms1_peaks);
  peak.idx <- 1;
  k <- 1;
  features.l <- list();
  
  while( peak.idx <= nrow(ms1_peaks.df) ){
    if( isFALSE(ms1_peaks.df[peak.idx, 'is_ion']) ){
      
      # print(peak.idx)
      
      ms1_mat <- extract_ms_matrix.f(peak.idx = peak.idx, ms1_peaks.df = ms1_peaks.df, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                     ppm.n = ppm.n, rt_index = TRUE, mz_range = NULL, iso_win_index = NULL);
      if( nrow(ms1_mat) <= 1 ){
        print(paste( peak.idx, 'No MS1 data.'))
        peak.idx <- peak.idx + 1
        next 
      };
      
      # NMF on MS1 data
      ms1_rank <- 3 # to be modified later on, it should be detected automatically
      ngmcas_res <- nGMCAs(X.m = ms1_mat, rank = ms1_rank,
                           maximumIteration = maximumIteration, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB,
                           initialization_method = MS1_init_method,
                           errors_print = errors_print);
      
      W_ms1 <- ngmcas_res$S
      H_ms1 <- ngmcas_res$A
      W_ms1 <- apply(W_ms1, 2, function(x) x / max(x));  # normalize every column (comp spectrum) by the max value
      
      # prepare pure MS1 spectra and choose the peak corresponding spectrum (the good spectrum)
      mz_prec <- as.numeric(ms1_peaks.df[peak.idx, 'mz']);
      rt_prec <- as.numeric(ms1_peaks.df[peak.idx, 'rt']);
      
      # choose the good MS1 spectrum
      ms1_pure_data <- extract_ms1_pure_spectrum(W_ms1 = W_ms1, mz_prec = mz_prec);
      ms1_pure_spectrum <- ms1_pure_data$ms1_pure_spectrum;
      comp_ms1 <- ms1_pure_data$comp_ms1;
      
      # Test the chosen ms1 pure spectra ions, which will also be considered as peaks or not.
      ions_are_peaks <- check_ms1_ions(W_ms1 = W_ms1, comp_ms1 = comp_ms1, ms1_peaks.df = ms1_peaks.df, rt_prec = rt_prec, rt_tol = rt_tol);
      # these ions will not factorized again, but they may be used in different peaks factorization
      ms1_peaks.df[ions_are_peaks, 'is_ion'] <- TRUE;
      # --------------------------------------------------------------------------------------------------------- the peaks data.frame is updated :).
      
      if ( ms_level == "MS1" ){
        
        feature_sub.l <- list( 
          'peak' = ms1_peaks.df[peak.idx, ],
          'MS1_pure_spectrum' = ms1_pure_spectrum,
          'MS1_mixed_mat' = ms1_mat,
          'W_ms1' = W_ms1,
          'H_ms1' = H_ms1,
          'comp_ms1' = comp_ms1,
          'ms1_ions_are_peaks' = ions_are_peaks  );
        
      } else { # the user want to process the MS2 data also
        
        info.swath <- isolationWindows.range(rawData.onDiskMSnExp);
        
        res_ms2 <- extract_ms2_matrices(peak.idx = peak.idx, ms1_peaks.df = ms1_peaks.df,
                                        ppm.n = ppm.n, ms1_pure_spectrum = ms1_pure_spectrum,
                                        rawData.onDiskMSnExp = rawData.onDiskMSnExp, info.swath = info.swath );
        if( length(res_ms2) == 0 ){
          print(paste( peak.idx,'No MS2 data.'))
          peak.idx <- peak.idx + 1
          next
        };
        ms2_matrices <- do.call(rbind, res_ms2);
        if( nrow(ms2_matrices) <= 1 ){
          print(paste( peak.idx, 'No enough MS2 data.'))
          peak.idx <- peak.idx + 1
          next
        };
  
        # NMF on MS2 data
        ms2_rank <- 3  # change this !!!!!!!!!
        rownames(H_ms1) <- NULL; 
        ngmcas_res_all <- nGMCAs(X.m = ms2_matrices, rank = ms2_rank,
                                 maximumIteration = maximumIteration, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB,
                                 initialization_method = MS2_init_method, H_sub = H_ms1,
                                 errors_print = errors_print);
        W_ms2 <- ngmcas_res_all$S;
        H_ms2 <- ngmcas_res_all$A;
        W_ms2 <- apply(W_ms2, 2, function(x) x / max(x));  # normalize every column (comp spectrum) by the max value
        
        # choose the good MS2 spectrum
        comp_ms2 <- elutions_corelation(chromo_main = H_ms1[comp_ms1, ], chromos = H_ms2);
        ms2_pure_spectrum <- choose_ms2_pure_spectrum(W_ms2 = W_ms2, choosen_comp = comp_ms2);
        
        # Filter the ms2 pure spectra to get only the fragments related to the precursor from the SWATH window where it was fragmented
        ms2_pure_spectrum_specific <- filter_ms2_spectrum(ms2_pure_spectrum = ms2_pure_spectrum, ms2_matrices = res_ms2, mz_prec, info.swath = info.swath, peak.idx);
        
        # now I can delete the zero intensity fragments
        ms2_pure_spectrum <- ms2_pure_spectrum[ms2_pure_spectrum['intensity'] != 0, ];
        
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
          'comp_ms1' = comp_ms1,
          'comp_ms2' = comp_ms2,
          'ions_are_peaks' = ions_are_peaks  );
        
      }
      
      # # the following line will be deleted at the end !!!!!!!!!!!!!!!!!!!!!!!!!!
      # saveRDS(feature_sub.l, paste0('C:/Users/DK273056/Documents/DIA_NMF_R_package_outputs/peak_', peak.idx, '.rds') );
      
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
