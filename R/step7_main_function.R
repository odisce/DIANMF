#' main function; input: some parameters and mzML file path, output: list of features.
#' 
#' @param mzML_path mzML file path.
#' @inheritParams extract_mixedMat
#' @param peaks_by_xcms `Logical` `TRUE` if the user wants to detect MS1 levels using XCMS, `FALSE` if the user will provide the peaks matrix or data frame to be identified.
#' @inheritParams prepare_ms1_peaks
#' @param d.out `character` file path to save the MS1 peaks matrix if detected by XCMS.
#' @param put_sn_thr `numeric` to put sn threshold on the peaks, else NULL.
#' @inheritParams nGMCAs
#' @param MS1_init_method `character` MS1 initialization method.
#' @param MS2_init_method `character` MS2 initialization method.
#' @inheritParams check_ms1_ions
#' @inheritParams find_rank
#' @inheritParams extract_ms_matrix.f
#' @inheritParams pure_sources.f
#' @param ... Additional parameters passed to \code{\link[xcms]{CentWaveParam}}.
#'
#' @return `list` of identified features.
#' 
#' @export
#' 
#' @import MSnbase
#' @import magrittr
#' @import dplyr 
#' @importFrom purrr map
#' @importFrom data.table %between%
dia_nmf.f <- function(
    mzML_path = NULL,
    ms_level = c("MS1" ,"MS2"),
    # parameters to detect peaks by xcms, or input the peak matrix or data.frame
    peaks_by_xcms = c(TRUE, FALSE), ms1_peaks = NULL, d.out = NULL, put_sn_thr = 3,
    # parameters to extract MS1 & MS2 mixed matrices
    ppm.n = 7, mz_range = NULL,
    # NMF parameters
    maximumIteration = 10, maxFBIteration = 10, toleranceFB = 1e-5,
    MS1_init_method = c('nndsvd', 'random'), MS2_init_method = c('nndsvd', 'subSample', 'random'), errors_print = FALSE,
    # rt tolerance 
    rt_tol_ions = 2, # used to find MS1 ions that also detect as MS1 peaks
    # pure sources parameters
    ms_type = 'sum',
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

  if( !is.null(put_sn_thr) ){  # Delete peaks of sn <= put_sn_thr
    ms1_peaks <- as.data.frame(ms1_peaks)
    ms1_peaks <- ms1_peaks[ms1_peaks$sn >= put_sn_thr, ]
    saveRDS(ms1_peaks, file = paste0(d.out, '/ms1_peaks.rds'))
  };

  ms1_peaks.df <- prepare_ms1_peaks(ms1_peaks = ms1_peaks);
  ms1_peaks.df$r <- seq(1, nrow(ms1_peaks.df));
  
  peak.idx <- 1;
  k <- 1;
  features.l <- list();
  
  while( peak.idx <= nrow(ms1_peaks.df) ){
    if( ms1_peaks.df[peak.idx, 'is_ion'] == 0 ){
      
      print(peak.idx)
      
      # check if it is a real peak or noise ------------------------------------
      peak <- ms1_peaks.df[peak.idx, ]
      rtr <- c(peak$rtmin-5, peak$rtmax+5)
      mzr <- c(peak$mzmin, peak$mzmax)
      xic_data <- rawData.onDiskMSnExp |>
        MSnbase::filterMsLevel(1L) |>
        MSnbase::filterRt(rt = rtr) |>
        MSnbase::filterMz(mz = mzr)
      
      rt <- unname(rtime(xic_data))      
      inten <- intensity(xic_data)
      for( x in 1:length(inten) ){
        if(length(inten[[x]]) == 0){
          inten[[x]] <- 0  }  }
      for( x in 1:length(inten) ){
        inten[[x]] <- sum(inten[[x]])  }
      
      inten <- unname(unlist(inten))
      
      chrom <- data.frame(
        'rt' = rt,
        'intensity' = inten )
      
      is.peak <- is_true_peak(chromatogram = chrom) 
      if( isFALSE(is.peak) ){
        print(paste( peak.idx, 'Not a real peak; noise.'))
        peak.idx <- peak.idx + 1
        next 
      };
      # ------------------------------------------------------------------------Done;
      
      mz_prec <- as.numeric(ms1_peaks.df[peak.idx, 'mz']);
      rt_prec <- as.numeric(ms1_peaks.df[peak.idx, 'rt']);

      ms1_mat <- extract_ms_matrix.f(peak.idx = peak.idx, ms1_peaks.df = ms1_peaks.df, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                     ppm.n = ppm.n, rt_index = TRUE, mz_range = mz_range, iso_win_index = NULL);
      if( is.null(ms1_mat) ){
        print(paste( peak.idx, 'No MS1 data.'))
        peak.idx <- peak.idx + 1
        next 
      };
      if( nrow(ms1_mat) <= 1 ){
        print(paste( peak.idx, 'No MS1 data.'))
        peak.idx <- peak.idx + 1
        next 
      };

      # determine the rank of factorization
      rt_axis <- as.numeric(colnames(ms1_mat))
      rank <-  find_rank(ms1_peaks.df, peak.idx, min_rt = rt_axis[[1]], max_rt = rt_axis[[length(rt_axis)]], max_r = ncol(ms1_mat));
      if( rank == 0 ){
        # print(paste(peak.idx, "No factorization"))
        peak.idx <- peak.idx + 1
        next
      };
      
      # NMF on MS1 data
      ngmcas_res <- nGMCAs(X.m = ms1_mat, rank = rank,
                           maximumIteration = maximumIteration, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB,
                           initialization_method = MS1_init_method,
                           errors_print = errors_print);
      
      if( is.null(ngmcas_res) ){ 
        peak.idx <- peak.idx + 1
        next
      };
      
      W_ms1 <- ngmcas_res$S;
      H_ms1 <- ngmcas_res$A;
      colnames(H_ms1) <- colnames(ms1_mat);
      
      # retrieve the MS1 sources
      ms1_pure_sources <- pure_sources.f(W = W_ms1, H = H_ms1, ms_type = ms_type);
      
      # prepare pure MS1 spectra and choose the peak corresponding spectrum (the good spectrum)
      # choose the good MS1 spectrum
      ms1_spectra_mat <- matrix(0, nrow = nrow(ms1_pure_sources[[1]]$source_spect), ncol = length(ms1_pure_sources));
      for(i in 1:length(ms1_pure_sources)){
        ms1_spectra_mat[, i] <- ms1_pure_sources[[i]]$source_spect$intensity
      };
      rownames(ms1_spectra_mat) <- ms1_pure_sources[[i]]$source_spect$mz_value;
      
      ms1_pure_data <- extract_ms1_pure_spectrum(W_ms1 = ms1_spectra_mat, mz_prec = mz_prec);
      ms1_pure_spectrum <- ms1_pure_data$ms1_pure_spectrum;
      comp_ms1 <- ms1_pure_data$comp_ms1;
      
      # Test the chosen ms1 pure spectra ions, which will also be considered as peaks or not.--------------------------------------------------------
      ions_are_peaks <- check_ms1_ions(W_ms1 = ms1_spectra_mat, comp_ms1 = comp_ms1, ms1_peaks.df = ms1_peaks.df, rt_prec = rt_prec, rt_tol_ions = rt_tol_ions);
      # these ions will not factorized again, but they may be used in different peaks factorization
      ms1_peaks.df[ions_are_peaks, 'is_ion'] <- peak.idx;
      
      ms1_peaks.df[peak.idx, "is_ion"] <- peak.idx; # no need, for sure the precursor is one of the ions peaks
      # --------------------------------------------------------------------------------------------------------- the peaks data.frame is updated :).
      
      if ( ms_level == "MS1" ){
        
        feature_sub.l <- list( 
          'peak' = ms1_peaks.df[peak.idx, ],
          'MS1_mixed_mat' = ms1_mat,
          # 'W_ms1' = W_ms1,
          # 'H_ms1' = H_ms1,
          'comp_ms1' = comp_ms1,
          'rank' = rank,
          'ms1_ions_are_peaks' = ions_are_peaks,
          'MS1_pure_spectrum' = ms1_pure_spectrum, 
          'ms1_pure_sources' = ms1_pure_sources  );
        
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
        rownames(H_ms1) <- NULL; 
        ngmcas_res_all <- nGMCAs(X.m = ms2_matrices, rank = rank,
                                 maximumIteration = maximumIteration, maxFBIteration = maxFBIteration, toleranceFB = toleranceFB,
                                 initialization_method = MS2_init_method, H_sub = H_ms1,
                                 errors_print = errors_print);
        if( is.null(ngmcas_res_all) ){ 
          peak.idx <- peak.idx + 1
          next
        };
        
        W_ms2 <- ngmcas_res_all$S;
        rownames(W_ms2) <- rownames(ms2_matrices);
        H_ms2 <- ngmcas_res_all$A;

        # retrieve the MS2 sources
        ms2_pure_sources <- pure_sources.f(W = W_ms2, H = H_ms2, ms_type = ms_type);
        
        # choose the good MS2 spectrum
        ms2_spectra_mat <- matrix(0, nrow = nrow(ms2_pure_sources[[1]]$source_spect), ncol = length(ms2_pure_sources));
        for(i in 1:length(ms2_pure_sources)){
          ms2_spectra_mat[, i] <- ms2_pure_sources[[i]]$source_spect$intensity
        };
        rownames(ms2_spectra_mat) <- ms2_pure_sources[[i]]$source_spect$mz_value;
        
        comp_ms2 <- elutions_corelation(chromo_main = H_ms1[comp_ms1, ], chromos = H_ms2);
        ms2_pure_spectrum <- choose_ms2_pure_spectrum(W_ms2 = ms2_spectra_mat, choosen_comp = comp_ms2);
        
        # Filter the ms2 pure spectra to get only the fragments related to the precursor from the SWATH window where it was fragmented
        # ms2_pure_spectrum_specific <- filter_ms2_spectrum(ms2_pure_spectrum = ms2_pure_spectrum, ms2_matrices = res_ms2, mz_prec, info.swath = info.swath);
        # ms2_pure_spectrum_specific$mz_value <- gsub("\\.\\d$", "", ms2_pure_spectrum_specific$mz_value); # Remove trailing '.X' where X is any digit
        # ms2_pure_spectrum_specific$mz_value <- as.numeric(ms2_pure_spectrum_specific$mz_value);
        
        # now I can delete the zero intensity fragments
        ms2_pure_spectrum <- ms2_pure_spectrum[ms2_pure_spectrum['intensity'] != 0, ];
        ms2_pure_spectrum$mz_value <- gsub("\\.\\d$", "", ms2_pure_spectrum$mz_value); # Remove trailing '.X' where X is any digit
        ms2_pure_spectrum$mz_value <- as.numeric(ms2_pure_spectrum$mz_value);
        
        feature_sub.l <- list(
          'peak' = ms1_peaks.df[peak.idx, ],
          'MS1_mixed_mat' = ms1_mat,
          'MS2_mixed_mat' = res_ms2,
          # 'W_ms1' = W_ms1,
          # 'H_ms1' = H_ms1,
          # 'W_ms2' = W_ms2,
          # 'H_ms2' = H_ms2,
          'comp_ms1' = comp_ms1,
          'comp_ms2' = comp_ms2,
          'ions_are_peaks' = ions_are_peaks,
          'rank' = rank,
          'ms1_pure_sources' = ms1_pure_sources,  # real sources
          'ms2_pure_sources' = ms2_pure_sources,
          'MS1_pure_spectrum' = ms1_pure_spectrum,   # real spectra
          'MS2_pure_spectrum' = ms2_pure_spectrum );
          # 'MS2_pure_spectrum_specific' = ms2_pure_spectrum_specific  );
        
      }

      features.l[[k]] <- feature_sub.l;
      k <- k + 1
      } else {   # peaks which are ions will not be factorized alone
          peak.idx <- peak.idx + 1
          next
        }
    
      peak.idx <- peak.idx + 1
  }
  
  if (!is.null(d.out)){
    saveRDS(ms1_peaks.df, file = paste0(d.out, '/ms1_peaks_updated.rds'))
  }
  
  return(features.l)
}
