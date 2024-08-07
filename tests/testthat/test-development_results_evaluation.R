skip("development script")

test_that("Evaluate the final results", {
  require(ggplot2)
  require(dplyr)
  require(parallel)
  NewDB <- readRDS("~/1workflow/DIA_NMF_workflow/Data/NewDB.rds")
  file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
  rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");
  # load the isolation windows
  iso <- isolationWindows.range(rawData.onDiskMSnExp)
  
  features <- readRDS("~/DIA_NMF_R_package_outputs/features.rds")
  peaks <- lapply(1:length(features), function(i){
    peak <- features[[i]]$peak
  })
  peaks <- do.call(rbind, peaks)
  
  # search for the spiked precursors in this peak list
  compounds.df <- load_compounds_function.tb("//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/Barbier_supplementary_tableS1.tsv")
  precursors_names <- readRDS("~/1workflow/DIA_NMF_workflow/Data/new_names.rds");
  precursors_names[5] <- "Triethanolamine ou Trolamine"
  precursors_names[34] <- "trans-Zeatin glucoside"
  rownames(compounds.df) <- precursors_names
  
  compounds.l <- create_list_function.l(compounds.df);
  l <- search_for_ROIs_function.l(compounds.l, peaks, rt_tol.n = 10)
  
  empty_elements <- sapply(l, function(x) length(x) == 0)
  l <- l[!empty_elements]
  prec_names <- names(l)
  
  # extract the pure ms2 spectra
  pure_spectra <- lapply(1:length(l), function(i){
    # i <- 5
    print(i)
    res <- features[[l[[i]][[1]]]]$MS2_pure_spectrum_specific
  })
  names(pure_spectra) <- names(l)
  
  # match the pure spectra
  scores <- lapply(1:length(pure_spectra), function(i){
    print(i)
    info_prec <- strsplit(prec_names[[i]], split = "_")
    mz_prec <- as.numeric(info_prec[[1]][[3]])
    x <- pure_spectra[[i]]
    sub_score <- match_pure_scores2(polarity = 'POS', mz_precursor = mz_prec, data_base = NewDB, measured_spectra = as.data.frame(x), mz_tol = 0.05)
    sub_score$true_name <- info_prec[[1]][[1]]
    row_idx <- which.max(sub_score$"total score");
    sub_score <-  sub_score[ row_idx, ]
  })
  scores <- do.call(rbind, scores)
  
  
  # plot them
  for (i in 1:length(l)) {
    lib_idx <- as.numeric(unname(scores[i, 'ref spectrum index']));
    lib_spect <- NewDB$spectra[[lib_idx]]@spectrum;
    colnames(lib_spect) <- c('mz_value', 'intensity');
    lib_spect <- as.data.frame(lib_spect);
    lib_spect$intensity <- lib_spect$intensity / max(lib_spect$intensity);
    my_spectrum <- pure_spectra[[i]]
    p1 <- ggplot2::ggplot( ) +
      geom_linerange(data = my_spectrum, stat = "identity", aes( x = mz_value, y = intensity, ymin = 0, ymax = intensity)) +
      geom_linerange(data = lib_spect, stat = "identity", aes( x = mz_value, y = -intensity, ymin = -intensity, ymax = 0, color = 'red')) +
      guides(color = FALSE) +
      labs( caption = paste('Dot prod=', scores[i,  "Dot prod"], ' Inverse prod=', scores[i,  "Rev prod"], ' Presence=', scores[i,  "Presence"], ';',  scores[i,  "name.name"], scores[i,  "ref spectrum index"])) +
      theme_bw();
    
    d.plot <- file.path('~/DIA_NMF_R_package_plots/pure_spectra');
    if (!dir.exists(d.plot)) {
      dir.create(d.plot, recursive = TRUE)
    };
    ggsave(paste0(d.plot, '/', prec_names[[i]], '.png'), p1, w=10, h = 10, dpi = 300)
  }
  
})


