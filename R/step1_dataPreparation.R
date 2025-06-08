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
#' @param mzml_dt `data.table` from `DIANMF::prepare_mzMLfiles()`
#' @return a data.table corresponding to a sequence.
#' 
#' @export
#' 
#' @import data.table
create_seq <- function(mzml_dt) {
  data.table(
    "mzml_path" = mzml_dt$mzml_path,
    "class" = "A",
    "InjectionOrder" = mzml_dt$InjectionOrder
  )
}


#' Detect MS1 peaks using XCMS.
#'
#' @param sequence_table A `data.frame` with at least (mzml_path, class, InjectionOrder).
#' @param params A `list` to set the parameters from xcms pipeline steps.
#'               List names should be xcms parameters methods (ex: `list("CentWaveParam" = CentWaveParam()`).
#' @return `MsExperiment` objects
#' @param rt_range (Optional) rt range to subset before running peak detection.
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
  output <- xcms::chromPeaks(msexp, isFilledColumn = TRUE) %>%
    data.table::as.data.table(keep.rownames = "peakid")
  output[, peakindex := seq_len(.N)]
  return(output[])
}


#' Extract MS1 features
#'
#' @inheritParams extract_xcms_peaks
#' @param orderL_sample `character` sample name.
#' @param quantifyL `Logical` to perform quantification with `xcms::quantify()`
#'
#' @return `data.table` `data.frame` object.
#' @export
#'
#' @importFrom xcms featureDefinitions quantify
#' @importFrom data.table as.data.table
#' @importFrom SummarizedExperiment assay
#' @importFrom MsExperiment sampleData
#' @import magrittr
extract_xcms_features <- function(msexp, orderL_sample = NULL, quantifyL = FALSE) {
  res <- MsExperiment::sampleData(msexp) %>% as.data.table()
  res <- basename(res$spectraOrigin)
  
  output <- msexp %>%
    xcms::featureDefinitions(.) %>%
    data.table::as.data.table(keep.rownames = "featureid")
  if (quantifyL) {
    intensities <- msexp %>%
      xcms::quantify(., method = "max") %>%
      SummarizedExperiment::assay() %>%
      as.data.table(., keep.rownames = "rn")
    colnames(intensities)[1] <- "featureid"
    output <- merge(output, intensities, by = "featureid")
    if (!is.null(orderL_sample)) {
      if (!orderL_sample %in% res) {
        stop(sprintf("sample not found: %s", orderL_sample))
      }
      output <- output[order(-get(orderL_sample)), ]
    }
  }
  return(output)
}