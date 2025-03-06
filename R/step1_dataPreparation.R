# Step1: load mzML files, detect and align features and peaks by xcms

#' Load and prepare mzML files
#'
#' @param input_dir `character` mzML files directory path.
#'
#' @return `data.table` `data.frame` contains mzML files information.
#' @export
#' 
#' @import data.table 
#' @importFrom tools file_path_sans_ext
prepare_mzMLfiles <- function(input_dir){
  input_files <- list.files(input_dir, pattern = ".mzml", full.names = TRUE, ignore.case = TRUE)
  mzml_dt <- data.table(
    "mzml_path" = input_files,
    "file_name" = tools::file_path_sans_ext(basename(input_files)),
    "creation_date" = file.info(input_files)$ctime %>% as.POSIXct(),
    "class" = gsub("^.*-([0-9]{1,3}(,[0-9]{1,3})?)ng.*$", "\\1", basename(input_files)) %>%
      gsub(",", ".", .) %>%
      as.numeric()
  )
  mzml_dt[order(creation_date), InjectionOrder := seq_len(.N)]
  
  return(mzml_dt)
}


#' Create simple sequence from a list of mzML paths
#'
#' @param mzML_path A character vector with mzML paths.
#' @return a data.table corresponding to a sequence.
#' 
#' @export
#' 
#' @import data.table
create_seq <- function(mzML_path) {
  data.table(
    "mzml_path" = mzML_path,
    "class" = "A",
    "InjectionOrder" = seq_len(length(mzML_path))
  )
}


#' Detect MS1 peaks using XCMS
#'
#' @param sequence_table A data.frame with at least (mzml_path, class, InjectionOrder).
#' @param params A list to set the parameters from xcms pipeline steps.
#'               List names should be xcms parameters methods (ex: `list("CentWaveParam" = CentWaveParam()`).
#' @param rt_range (Optional) retention time range to subset before running peak detection.
#'
#' @return `MsExperiment` objects.
#' 
#' @export
#' 
#' @import MsExperiment MsFeatures xcms data.table magrittr
detect_xcms_peaks <- function(
    sequence_table,
    params = list(
      "CentWaveParam" = xcms::CentWaveParam(),
      "MergeNeighboringPeaksParam" = xcms::MergeNeighboringPeaksParam(),
      "ObiwarpParam" = xcms::ObiwarpParam(),
      "PeakDensityParam" = xcms::PeakDensityParam(sampleGroups = NA),
      "ChromPeakAreaParam" = xcms::ChromPeakAreaParam()
    ),
    rt_range = NULL
){
  sequence_table <- data.table::as.data.table(sequence_table)
  sequence_table <- sequence_table[order(InjectionOrder), ]
  xcms_obj <- MsExperiment::readMsExperiment(
    spectraFiles = sequence_table$mzml_path,
    sampleData = sequence_table
  )
  if (!is.null(rt_range)) {
    xcms_obj <- xcms::filterRt(xcms_obj, rt_range)
  }
  if (length(params$PeakDensityParam@sampleGroups) == 1 && is.na(params$PeakDensityParam@sampleGroups)) {
    params$PeakDensityParam@sampleGroups <- sampleData(xcms_obj)$class
  }
  xcms_obj_peaks <- xcms::findChromPeaks(
    xcms_obj,
    param = params$CentWaveParam
  ) %>% 
    refineChromPeaks(., param = params$MergeNeighboringPeaksParam, msLevel = 1L)
  if (length(xcms::fileNames(xcms_obj_peaks)) > 1) {
    xcms_obj_peaks <- adjustRtime(xcms_obj_peaks, param = params$ObiwarpParam)
  }
  output <- groupChromPeaks(xcms_obj_peaks, param = params$PeakDensityParam)
  
  if (length(xcms::fileNames(xcms_obj_peaks)) > 1) {
    output <- fillChromPeaks(output, param = params$ChromPeakAreaParam)
  }
  return(output)
}


#' Extract MS1 peaks
#'
#' @param msexp `MsExperiment` object obtained from xcms or with `DIANMF::detect_xcms_peaks()`.
#'
#' @return MS1 peaks `matrix`.
#' 
#' @export
#' 
#' @importFrom xcms findChromPeaks CentWaveParam chromPeaks
#' @import magrittr MsExperiment data.table
extract_xcms_peaks <- function(msexp) {
  output <- xcms::chromPeaks(msexp) %>%
    data.table::as.data.table(keep.rownames = "peakid")
  return(output)
}


#' Detect LC features
#'
#' @param mzml_dt `data.table` created in` DIANMF::prepare_mzMLfiles`
#' @param temp_saveL `logical`
#' @inheritParams detect_xcms_peaks
#'
#' @return `XcmsExperiment` `xcms` object.
#' @export
detect_LCfeatures <- function(mzml_dt, params, temp_saveL = T, rt_range){
  
  if (temp_saveL) {
    save_path <- "./temp/data/"
    dir.create(save_path, recursive = TRUE)
    xcms_obj_path <- file.path(save_path, "xcms_obj.rds")
    if (file.exists(xcms_obj_path)) {
      xcms_obj <- readRDS(xcms_obj_path)
    } else {
      xcms_obj <- detect_xcms_peaks(
        sequence_table = mzml_dt,
        params,
        rt_range
      )
      saveRDS(xcms_obj, file = xcms_obj_path)
    }
  }
  
  return(xcms_obj)
}


#' Extract MS1 features
#'
#' @inheritParams extract_xcms_peaks
#' @param orderL `logical` TRUE to order the features by intensity order of orderL_sample.
#' @param orderL_sample `character` sample name.
#'
#' @return `data.table` `data.frame` object.
#' @export
#'
#' @importFrom xcms featureDefinitions quantify
#' @importFrom data.table as.data.table
#' @importFrom SummarizedExperiment assay
#' @import magrittr
extract_xcms_features <- function(msexp, orderL = F, orderL_sample = NULL) {
  
  output <- msexp %>%
    xcms::featureDefinitions(.) %>%
    data.table::as.data.table(keep.rownames = "featureid")
  
  intensities <- msexp %>% 
    xcms::quantify(., method = "max") %>%
    SummarizedExperiment::assay() %>%
    as.data.table(., keep.rownames="rn")
  colnames(intensities)[1] <- "featureid"
  
  res <- merge(output, intensities, by = "featureid")
  
  if (orderL & !is.null(orderL_sample)) {
    res <- res[order(-get(orderL_sample)), ]
  }
  
  return(res)
}

