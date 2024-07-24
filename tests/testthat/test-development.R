skip("development script")

test_that("test on barbier data 10ng replicate 3", {
  require(xcms)
  require(MSnbase)
  require(ggplot2)
  require(dplyr)
  require(parallel)
  
  devtools::load_all()
  
  file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
  rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");
  fast_sample <- FALSE
  
  if (fast_sample) {
    sample_file <- MSnbase::filterRt(rawData.onDiskMSnExp, c(0,2*60))
  } else {
    sample_file <- rawData.onDiskMSnExp
  }
  eics_peaks <- detect_EICs.f(
    rawData.onDiskMSnExp = sample_file,
    d.out = NULL,
    ppm = 6,
    peakwidth = c(6,60), 
    snthresh = 1,
    prefilter = c(5,4000), 
    mzCenterFun = "wMeanApex3",
    integrate = 2,
    mzdiff = -0.001, 
    noise = 0, 
    firstBaselineCheck = FALSE
  )
  eics_peaks.mat <- xcms::chromPeaks(eics_peaks)
  
  expect_true("XCMSnExp" %in% class(eics_peaks))
  expect_true("matrix" %in% class(eics_peaks.mat))
  expect_true( nrow(eics_peaks.mat) >= 1 )
  
  eics_peaks.mat <- as.matrix(eics_peaks.mat[order(-eics_peaks.mat[, 'into']), ])
  idx.pg <- 1

  ms1_mat <- extract_ms_matrix.f(idx.pg = idx.pg, eics_peaks.mat = eics_peaks.mat, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                  ppm = 7, rt_index = TRUE, mz_range = NULL, iso_win_index = NULL)
   
  expect_true("matrix" %in% class(ms1_mat))
  expect_true(dim(ms1_mat)[1] & dim(ms1_mat)[2] > 0)

  # ms2_mat <- 
  
  # expect_true("matrix" %in% class(ms2_mat))
  # expect_true(dim(ms2_mat)[1] & dim(ms2_mat)[2] > 0)
  # expect_equal(ncol(ms1_mat), ncol(ms2_mat))
  
  # plot raw data: before factorization
  # plot the mixed MS1 eics
  mat <- ms1_mat
  ms1_rt <- as.numeric(colnames(mat))
  ms1_mz_values <- paste0("MS1_", as.numeric(rownames(mat)))
  mat <- matrix(mat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(mat) <- ms1_mz_values
  eics <- prepare_mixed_eics(eics = mat, retention_times = ms1_rt)
  ggplot2::ggplot( data = eics, aes(x = rt, y = intensity, color = MS_type)) +
    geom_vline( xintercept = as.numeric( 690.44971 ), colour = "black", linetype = 1, size = 0.5 ) +
    geom_line() +
    geom_point() +
    # geom_line(data = subset(eics, MS_type == "MS1"), color = "black", size = 0.8, linetype = 2) +
    xlim( min(eics$rt), max(eics$rt) ) +
    guides(color = FALSE)
  
  # plot the mixed MS1 spectrum
  ms1_mixed_spectrum <- prepare_mixed_spectra(mat, as.numeric(rownames(ms1_mat)))
  ms1_mixed_spectrum$intensity <- ms1_mixed_spectrum$intensity / max(ms1_mixed_spectrum$intensity)
  ggplot2::ggplot( ) +
    geom_linerange(data = ms1_mixed_spectrum, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
    # guides(color = FALSE) +
    labs( caption = paste("MS1 mixed spectrum of the most intense peak of mz", round(mz_prec,4))) +
    theme_bw();
  
  
  # Factorize ms1_mat using nGMCAs
  options.l <-  list(
    'maximumIteration' = 10,
    'maxFBIteration' = 5,
    'toleranceFB' = 1e-5,
    'useTranspose' = TRUE)
  ngmcas_res <- nGMCAs(originalMatrix = ms1_mat, rank = 2, options = options.l, errors_print = TRUE, initialization_method = 'nndsvd')
  W <- ngmcas_res$S
  H <- ngmcas_res$A

  # extract the corresponding pure spectra: prepare and plot
  mz_prec <- eics_peaks.mat[idx.pg, 'mz']
  rt_prec <- eics_peaks.mat[idx.pg, 'rt']
  
  mz_ms1_ions <- round(as.numeric(rownames(W)), 4);
  rownames(W) <- NULL
  spectra <- prepare_pure_spectra(W , mz_ms1_ions)
  
  # choose the component
  closest_row <- which.min(abs(mz_ms1_ions - mz_prec));
  choosen_comp <- which.max(W[closest_row, ]);
  spectrum <- spectra[spectra$comp_nb == paste0('comp',choosen_comp), ]
  spectrum$intensity <- spectrum$intensity / max(spectrum$intensity)
  
  # plot the MS1 pure spectrum of this peak
  ggplot2::ggplot( ) +
    geom_linerange(data = spectrum, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
    # guides(color = FALSE) +
    labs( caption = paste("MS1 pure spectrum of the most intense peak of mz", round(mz_prec,4))) +
    theme_bw();
  
  # this pure MS1 spectra supposed to contains all MS1 ions related to this peak; i.e. adducts, losses, isotopes and fragments
  # I will check fragments of min and max mz in which isolation windows belong
  min_mz <- min(spectrum$mz_value)
  max_mz <- max(spectrum$mz_value)
  info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
  idx.swath <- which(info.swath$lowerMz %between% c(min_mz, max_mz) | info.swath$upperMz %between% c(min_mz, max_mz)) # check this for other values
    
  # so from every isolation window I will extract the mixed ms2 data at the apex of the MS1 peak
  res <- lapply(idx.swath, function(i){
    # i <- 1
    ms2_mat <- extract_ms_matrix.f(idx.pg = idx.pg, eics_peaks.mat = eics_peaks.mat, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                   ppm = 7, rt_index = TRUE, mz_range = NULL, iso_win_index = i)
   
    return(ms2_mat)
  })
  # sapply(res, image) 
  # #------------------------------------------------------------------------------
  # require(gridExtra)
  # 
  # plist <- list()
  # 
  # for (i in 1:length(res)) {
  #   # i <- 1
  #   x <- res[[i]]
  #   x_rt <- as.numeric(colnames(x))
  #   x_mz_values <- paste0("MS2_", as.numeric(rownames(x)))
  #   eics <- matrix(x, nrow = nrow(x), ncol = ncol(x))
  #   rownames(eics) <- x_mz_values
  #   eics <- prepare_mixed_eics(eics = eics, retention_times = x_rt)
  #   
  #   plist[[i]] <- ggplot2::ggplot( data = eics, aes(x = rt, y = intensity, color = MS_type)) +
  #     geom_vline( xintercept = as.numeric( 690.44971 ), colour = "black", linetype = 1, size = 0.5 ) +
  #     geom_line() +
  #     geom_point() +
  #     # geom_line(data = subset(eics, MS_type == "MS1"), color = "black", size = 0.8, linetype = 2) +
  #     xlim( min(eics$rt), max(eics$rt) ) +
  #     guides(color = FALSE)
  #   
  #   pdf(sprintf("p%s.pdf", i),
  #       width = 10, height = 4, onefile = T)
  #   plot(plist[[i]])
  #   dev.off()
  # }
  # 
  # do.call(grid.arrange, c(plist, ncol = 10))
  # # ---------------------------------------------------------------------------------------
  # plist <- list()
  # y_min <- Inf
  # y_max <- -Inf
  # 
  # # First loop to determine the y-axis limits
  # for (i in 1:length(res)) {
  #   x <- res[[i]]
  #   x_rt <- as.numeric(colnames(x))
  #   x_mz_values <- paste0("MS2_", as.numeric(rownames(x)))
  #   eics <- matrix(x, nrow = nrow(x), ncol = ncol(x))
  #   rownames(eics) <- x_mz_values
  #   eics <- prepare_mixed_eics(eics = eics, retention_times = x_rt)
  #   y_min <- min(y_min, min(eics$intensity))
  #   y_max <- max(y_max, max(eics$intensity))
  # }
  # 
  # # Second loop to create the plots with consistent y-axis limits
  # for (i in 1:length(res)) {
  #   x <- res[[i]]
  #   x_rt <- as.numeric(colnames(x))
  #   x_mz_values <- paste0("MS2_", as.numeric(rownames(x)))
  #   eics <- matrix(x, nrow = nrow(x), ncol = ncol(x))
  #   rownames(eics) <- x_mz_values
  #   eics <- prepare_mixed_eics(eics = eics, retention_times = x_rt)
  #   
  #   p <- ggplot2::ggplot(data = eics, aes(x = rt, y = intensity, color = MS_type)) +
  #     geom_vline(xintercept = as.numeric(690.44971), colour = "black", linetype = 1, size = 0.5) +
  #     geom_line() +
  #     geom_point() +
  #     xlim(min(eics$rt), max(eics$rt)) +
  #     ylim(y_min, y_max) +  # Set consistent y-axis limits
  #     guides(color = FALSE)
  #   
  #   if (i != 1) {
  #     p <- p + theme(axis.title.y = element_blank(), 
  #                    axis.text.y = element_blank(), 
  #                    axis.ticks.y = element_blank())
  #   }
  #   
  #   plist[[i]] <- p
  #   
  #   pdf(sprintf("p%s.pdf", i), width = 10, height = 4, onefile = T)
  #   plot(plist[[i]])
  #   dev.off()
  # }
  # do.call(grid.arrange, c(plist, ncol = 10))
  # #-----------------------------------------------------------------------------
  # # Combine all datasets into one data frame
  # combined_eics <- data.frame()
  # for (i in 1:length(res)) {
  #   x <- res[[i]]
  #   x_rt <- as.numeric(colnames(x))
  #   x_mz_values <- paste0("MS2_", as.numeric(rownames(x)))
  #   eics <- matrix(x, nrow = nrow(x), ncol = ncol(x))
  #   rownames(eics) <- x_mz_values
  #   eics <- prepare_mixed_eics(eics = eics, retention_times = x_rt)
  #   eics$dataset <- paste0("dataset_", i)  # Add an identifier for each dataset
  #   combined_eics <- bind_rows(combined_eics, eics)
  # }
  # 
  # # Create the ggplot
  # p <- ggplot(combined_eics, aes(x = rt, y = intensity, color = dataset)) +
  #   geom_vline(xintercept = as.numeric(690.44971), colour = "black", linetype = 1, size = 0.5) +
  #   geom_line() +
  #   geom_point() +
  #   xlim(min(combined_eics$rt), max(combined_eics$rt)) +
  #   theme(legend.title = element_blank())
  # 
  # # Save the plot to a PDF
  # pdf("combined_plot.pdf", width = 10, height = 4)
  # print(p)
  # dev.off()
  # #-----------------------------------------------------------------------------
  
  # concatenate MS1 and MS2 matrices & apply NMF
  ms2_matrices <- do.call(rbind, res)
  all_matrices <- rbind(ms1_mat, ms2_matrices)
  
  # I should change the algorithm, add the initialization method!!!!!!!!!!!!!!!!!!
  ngmcas_res_all <- nGMCAs(originalMatrix = all_matrices, rank = 4, options = options.l, errors_print = TRUE, initialization_method = 'nndsvd')
  W <- ngmcas_res_all$S
  H <- ngmcas_res_all$A
  
  # # plot the pure elution profiles
  # pure_eics <- prepare_pure_eics(H, as.numeric(colnames(res[[1]])))
  # ggplot2::ggplot(data = pure_eics, aes(rt, intensity , color = comp_nb)) +
  #         facet_grid(comp_nb~.) +
  #         geom_line() +
  #         geom_point() +
  #         xlim(min(pure_eics$rt), max(pure_eics$rt))

  ###
  ms1_mz <- as.numeric(rownames(ms1_mat));
  ms2_mz <- as.numeric(rownames(ms2_matrices))
  mz_values <- c(ms1_mz, ms2_mz);
  
  numeric_values <- round(as.numeric(rownames(W)), 4);
  closest_row <- which.min(abs(numeric_values - mz_prec));
  choosen_comp <- which.max(W[closest_row, ]);
  
  rownames(W) <- NULL
  pure_spectra <- prepare_pure_spectra(W, mz_values)
  
  pure_spectra <-  pure_spectra %>%
    group_by(comp_nb) %>%
    slice(-(1:nrow(ms1_mat))) %>%
    mutate(intensity = intensity / max(intensity)) 
  
  choosen_spectra <- subset(pure_spectra, comp_nb == paste0('comp',choosen_comp));
  choosen_spectra <- choosen_spectra[, c(3,4)]
  # # plot the pure spectra
  # ggplot2::ggplot( ) +
  #   geom_linerange(data = choosen_spectra, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
  #   # guides(color = FALSE) +
  #   labs( caption = paste("MS2 pure spectrum related to this MS1 peak", round(mz_prec,4))) +
  #   theme_bw();
  # 
  # # plot all the pure spectra
  # ggplot2::ggplot( ) +
  #   geom_linerange(data = pure_spectra, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
  #   facet_grid(comp_nb~.) +
  #   # guides(color = FALSE) +
  #   labs( caption = paste("MS2 pure spectra", round(mz_prec,4))) +
  #   theme_bw();
  
  # update the peaks list
  # check which MS1 ions in this pure MS1 spectrum are also incliuded as peaks in the peak list
  ions_peak <- check_ms1_ions(spectrum.df = spectrum, peaks.mat = eics_peaks.mat, rt_prec = rt_prec, rt_tol = 0.1)
  ions_peaks_nb <- length(ions_peak)  # nb of MS1 ions which are MS1 peaks in the peak matrix/list
  # these ions will be deleted from the list of peaks
  # update the MS1 peaks.mat
  eics_peaks.mat <- eics_peaks.mat[-ions_peak, ]
  
  # save the MS1 & MS2 information in the feature list for this peak and go for the second peak
  
  
  # match the pure MS2 spectra with the data base
  # my_spectrum <- as.data.frame(choosen_spectra);
  # # delete the fragments of intensity 0
  # my_spectrum <- my_spectrum[my_spectrum$intensity != 0, ]
  # # if( max(my_spectrum$intensity) == 0 )  { return(NULL) }
  # 
  # NewDB <- readRDS("~/1workflow/DIA_NMF_workflow/Data/NewDB.rds")
  # scores.df <- match_pure_scores2(polarity = 'POS', mz_precursor = mz_prec, data_base = NewDB, measured_spectra = my_spectrum, mz_tol = NULL)
  # if( is.null(scores.df) ) { return(NULL) };
  # 
  # row_idx <- which.max(scores.df$"total score");
  # lib_idx <- as.numeric(unname(scores.df[row_idx, 'ref spectrum index']));
  # lib_spect <- NewDB$spectra[[lib_idx]]@spectrum;
  # colnames(lib_spect) <- c('mz_value', 'intensity');
  # lib_spect <- as.data.frame(lib_spect);
  # lib_spect$intensity <- lib_spect$intensity / max(lib_spect$intensity);
  # dot_prod <- as.numeric(unname(scores.df[row_idx, 'Dot prod']));
  # inverse_prod <- as.numeric(unname(scores.df[row_idx, 'Rev prod']));
  # prese <- as.numeric(unname(scores.df[row_idx, 'Presence']));
  # 
  # if( dot_prod == 0 & inverse_prod == 0 & inverse_prod == 0 ){
  #   scores.df <- NULL
  #   return(NULL) }
  # 
  # name_matched <- scores.df[row_idx, 'name.name'];
  # ggplot2::ggplot( ) +
  #   geom_linerange(data = my_spectrum, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
  #   geom_linerange(data = lib_spect, stat = "identity", aes( x = mz_value, y = -intensity, ymin = -intensity, ymax = 0, color = 'red')) +
  #   guides(color = FALSE) +
  #   labs( caption = paste('Dot prod=', dot_prod, ' Inverse prod=', inverse_prod, ' Presence=', prese, ';', name_matched )) +
  #   theme_bw();
  # # theme(text = element_text(size = 8));
  
  
  
})
