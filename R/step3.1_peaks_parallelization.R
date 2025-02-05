# file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
# rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");
# 
# eics_peaks.mat <- extract_xcms_peaks(rawData.onDiskMSnExp = rawData.onDiskMSnExp,
#                                        ppm = 7, peakwidth = c(6,60), snthresh = 1,
#                                        prefilter = c(5,4000), mzCenterFun = "wMeanApex3",
#                                        integrate = 2, mzdiff = -0.001, noise = 0,
#                                        firstBaselineCheck = FALSE )
# eics_peaks.df <- as.data.frame(eics_peaks.mat)
# 
# df_binned <- eics_peaks.df %>%
#   dplyr::mutate(rt_minute = floor(rt / 60)) %>%  # group rt by minutes
#   group_by(rt_minute)
# 
# rt_min <- unique(df_binned$rt_minute)
# res <- lapply(0:(length(rt_min)-1), function(i){
#   sub_df <- df_binned[df_binned$rt_minute == i, ]
# })
# 
# # total_res2 <- dia_nmf.f(mzML_path = file,
# #                         ms_level = "MS2",
# #                         peaks_by_xcms = FALSE, ms1_peaks = eics_peaks.mat)
# 
# 
# library(future.apply)
# plan(multisession, workers = 3)
# # plan(sequential) 
# processing_res <- future_lapply(1:length(res), function(i){
# 
#   proc_sub_res <- dia_nmf.f(mzML_path = file,
#                             ms_level = "MS2",
#                             peaks_by_xcms = FALSE, ms1_peaks = as.matrix(res[[i]]) )
# })
# 
# 
# processing_res <- future_lapply(1:length(res), function(i){
#   proc_sub_res <- dia_nmf.f(mzML_path = file,
#                             ms_level = "MS2",
#                             peaks_by_xcms = FALSE, 
#                             ms1_peaks = as.matrix(res[[i]]) )
#   saveRDS(proc_sub_res, paste0('C:/Users/DK273056/Documents/DIA_NMF_R_package_outputs2/', i, '.rds'))
# })


# processing_res <- future_lapply(1:length(res), function(i){
#   tryCatch({
#     proc_sub_res <- dia_nmf.f(mzML_path = file,
#                               ms_level = "MS2",
#                               peaks_by_xcms = FALSE, 
#                               ms1_peaks = as.matrix(res[[i]]) )
#   }, error = function(e) {
#     message(sprintf("Error in iteration ", i, e$message))
#     return(NULL)
#   })
# })
# not working!!!!!!!!!!!!!!!!!!

#--------------------------------------------------------------------------------
# library(datasets)
# library(stats)
# y <- lapply(mtcars, FUN = mean, trim = 0.10)
# 
# library(future.apply)
# plan(multisession) ## Run in parallel on local computer
# library(datasets)
# library(stats)
# y <- future_lapply(mtcars, FUN = mean, trim = 0.10)
# 
# microbenchmark::microbenchmark(
#   
#   "lapply" = {
#     lapply(mtcars, FUN = mean, trim = 0.10)
#     },
#   "future_lapply" = {
#     future_lapply(mtcars, FUN = mean, trim = 0.10)
#   },
#   times = 50,
#   unit = "seconds"
# )

# it is not faster!!! why?!

