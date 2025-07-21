# the following function, extract for every identified pure source in all iterations one Spectra object;
# this object contains 1 MS1 pure spectra and the 10 MS2 pure spectra from the different isolation windows.

#' Export MS Spectra
#'
#' @inheritParams get_elutionprofile
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table magrittr
get_feature_summary <- function(
  features.l,
  max_method = c("contribution", "max_value")[1]
) {
  feature_dt <- data.table()
  for (spli in seq_len(length(features.l))) {
    for (iteri in seq_len(length(features.l[[spli]]$PureFeatures))) {
      if (is.null(features.l[[spli]]$PureFeatures[[iteri]])) {
        next
      }
      out_dt <- features.l[[spli]]$PureFeatures[[iteri]]$ms_info[!is.na(source),]
      out_dt[, iteration := as.integer(iteri)]
      feature_dt <- rbind(feature_dt, out_dt)
      rm(out_dt)
    }
  }
  if (max_method == "max_value") {
    output <- feature_dt[order(-apex_val), .SD[which.max(apex_val), ], by = .(featureid, sample)]
  } else if (max_method == "contribution") {
    output <- feature_dt[order(-apex_val), .SD[which.max(contribution), ], by = .(featureid, sample)]
  } else {
    stop(sprintf("max_method not recognized: %s", max_method))
  }
  
  return(output[])
}


#' Export MS Spectra
#'
#' @inheritParams get_elutionprofile
#'
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table
get_feature_coord <- function(
  features.l,
  summary_dt = NULL,
  feature_id,
  sample_index = NULL,
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l, max_method = max_method)
  }
  feat_coord <- summary_dt[featureid == feature_id, ]
  if (!is.null(sample_index)) {
    feat_coord <- feat_coord[sample == sample_index, ]
  }
  if (nrow(feat_coord) == 0) {
    message(sprintf("Feature %s not extracted in sample: %i", feature_id, sample_index))
    return(FALSE)
  }
  return(feat_coord)
}

#' Extract elution profile for a specific feature
#'
#' @param features.l `list` obtained with `DIANMF::DIANMF.f()`
#' @param summary_dt `data.table` of features with their optimal sources obtained from DIANMF::get_feature_summary().
#' @param feature_id `character` feature id.
#' @param sample_index `numeric(1)` sample index.
#' @param type `string` to retrieve the "pure" or "mixed" elution profiles.
#' @param method `string` to extract "all" sources or just the "best" source.
#' @param max_method `string` to choose a method to select the optimal iteration:
#'   - `contribution`: select the first iteration with maximum contribution.
#'   - `max_value`: select the first iteration with maximum apex pure value.
#'
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table magrittr
get_elutionprofile <- function(
  features.l = features,
  summary_dt = NULL,
  feature_id,
  sample_index = NULL,
  type = c("pure", "mixed")[1],
  method = c("all", "best")[2],
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l, max_method = max_method)
  }
  feat_coord <- get_feature_coord(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    max_method = max_method
  )
  if (isFALSE(feat_coord)) {
    return(FALSE)
  }

  ## Extract MS1 / MS2
  if (type == "mixed") {
    mixed_names <- c("MS1_mixed_mat", "MS2_mixed_mat")
    rt_names <- c("MS1_pure_elution_profiles", "MS2_pure_elution_profiles")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i], rt_names[i]) %in% i_names)) {
        mixed_mat_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]] %>%
          melt() %>%
          as.data.table()
        rt_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[rt_names[i]]][, .(scan_norm, rtime)] %>%
          unique()
        setnames(mixed_mat_dt, c("xic_label", "scan_norm", "value"))
        output <- merge(mixed_mat_dt, rt_dt, by = "scan_norm")
        output[, mslevel := mslve_vc[i]]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
    mixed_mat_out <- mixed_mat_out[, .(xic_label, rtime, value, mslevel, rank = as.integer(NA))]
  } else if (type == "pure") {
    mixed_names <- c("MS1_pure_elution_profiles", "MS2_pure_elution_profiles")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i]) %in% i_names)) {
        source_to_get <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][, unique(rank)]
        if (method == "best") {
          source_to_get <- feat_coord[, source]
        }
        output <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][rank %in% source_to_get,]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
    mixed_mat_out <- mixed_mat_out[, .("xic_label" = rank, rtime, value, "mslevel" = MSid, rank = rank)]
  } else {
    stop('type argument not recognized')
  }
  return(mixed_mat_out)
}

#' Extract Spectra
#'
#' @inheritParams get_elutionprofile
#'
#' @returns `data.table` of features with their optimal sources
#' @export
#'
#' @import data.table magrittr
get_spectra <- function(
  features.l,
  summary_dt = NULL,
  feature_id,
  sample_index = 1,
  type = c("pure", "mixed")[1],
  method = c("all", "best")[2],
  max_method
) {
  if (is.null(summary_dt)) {
    summary_dt <- get_feature_summary(features.l = features.l, max_method = max_method)
  }
  feat_coord <- get_feature_coord(
    features.l = features.l,
    summary_dt = summary_dt,
    feature_id = feature_id,
    sample_index = sample_index,
    max_method = max_method
  )
  if (isFALSE(feat_coord)) {
    return(FALSE)
  }
  ## Extract MS1 / MS2
  if (type == "mixed") {
    mixed_names <- c("MS1_mixed_mat", "MS2_mixed_mat")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i], "ms_info") %in% i_names)) {
        mixed_mat_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]] %>%
          melt() %>%
          as.data.table()
        setnames(mixed_mat_dt, c("xic_label", "scan_norm", "value"))
        mixed_mat_dt <- mixed_mat_dt[, .(value = max(value)), by = .(xic_label)]
        ms_info_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][["ms_info"]]
        output <- merge(
          mixed_mat_dt,
          ms_info_dt[, .(xic_label, msLevel, IsoWin, featureid, mz, rt, rank = as.integer(NA))],
          by = "xic_label"
        )
        output[, mslevel := mslve_vc[i]]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
  } else if (type == "pure") {
    mixed_names <- c("ms1_pure_spectra", "ms2_pure_spectra")
    mslve_vc <- c("MS1", "MS2")
    mixed_mat_out <- data.table()
    for (i in seq_len(2)) {
      i_names <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]] %>% names()
      if (all(c(mixed_names[i]) %in% i_names)) {
        source_to_get <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][, unique(rank)]
        if (method == "best") {
          source_to_get <- feat_coord[, source]
        }
        output <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][[mixed_names[i]]][rank %in% source_to_get,]
        ms_info_dt <- features.l[[feat_coord$sample]]$PureFeatures[[feat_coord$iteration]][["ms_info"]]
        output <- merge(
          output,
          ms_info_dt[, .(xic_label, msLevel, IsoWin, featureid, mz, rt)],
          by = "xic_label"
        )
        output[, mslevel := mslve_vc[i]]
        mixed_mat_out <- rbind(mixed_mat_out, output, fill = TRUE)
        rm(output)
      }
    }
  } else {
    stop('type argument not recognized')
  }
  
  return(mixed_mat_out[, .(xic_label, IsoWin, mslevel, mz, value, rank)])
}

#' Extract spectrum (of class `Spectra`) from `data.table` `data.frame`
#'
#' @param ms_Spectrum `data.table` `data.frame`
#' @inheritParams get_spectra
#'
#' @returns `Spectra` object
#' @importFrom Spectra Spectra
#' @importFrom S4Vectors DataFrame
get_spectrum_obj <- function(ms_Spectrum, feature_id){
  
  # MS1
  ms1_df <- ms_Spectrum[ms_Spectrum$mslevel == "MS1", ]
  ms1_spd <- DataFrame(
    msLevel = 1L,
    name = feature_id,
    id = unique(ms1_df$rank),
    COMMENT = paste(ms1_df$xic_label, collapse = ";")  )
  
  ms1_spd$mz <- list(ms1_df$mz)
  ms1_spd$intensity <- list(ms1_df$value)
  ms1_spectrum <- Spectra(ms1_spd)
  
  # MS2
  ms2_df <- ms_Spectrum[ms_Spectrum$mslevel == "MS2", ]
  if( nrow(ms2_df) > 0 ){
      isolation_spectra <- lapply(unique(ms2_df$IsoWin), function(targetMz) {
      subset_df <- ms2_df[ms2_df$IsoWin == targetMz, ]  
      
      ms2_spd <- S4Vectors::DataFrame(
        msLevel = 2L,
        name = feature_id,
        id = unique(ms2_df$rank),
        isolationWindowTargetMz = as.numeric(targetMz)
        # isolationWindowLowerMz = as.numeric(unique(subset_df$isolationWindowLowerMz)),
        # isolationWindowUpperMz = as.numeric(unique(subset_df$isolationWindowUpperMz))
      )
      
      ms2_spd$mz <- list(subset_df$mz)
      ms2_spd$intensity <- list(subset_df$value)
      return(Spectra::Spectra(ms2_spd))
    })
    isolationWin_spectra <- do.call(c, isolation_spectra)
  } else {
    isolationWin_spectra <- NULL
  }
  
  ms_spectra <- c(ms1_spectrum, isolationWin_spectra)
  return(ms_spectra)
}

#' Extract `Spectra` object of a feature
#'
#' @inheritParams plot_Spectra
#'
#' @returns `Spectra` object of a feature
#' @export
export_featureSpect <- function(features.l, feature_id, sample_index, type = c("pure", "mixed")[1], method = c("all", "best")[2], max_method = c("contribution", "max_value")[2]){
  
  temp_ft <- get_feature_summary(features.l, max_method)
  temp_ft_sample <- temp_ft[ sample == sample_index, ]
  ms_Spectrum <- get_spectra(features.l, summary_dt = temp_ft_sample, feature_id,
                             sample_index,
                             type,
                             method,
                             max_method)
  
  if( nrow(ms_Spectrum) > 0){
    if( method == "best" ){
      sp <- get_spectrum_obj(ms_Spectrum, feature_id)
      return(sp)
    }else if(method == "all"){
      sp_diff_sources <- lapply(unique(ms_Spectrum$rank), function(source){
        sub_ms_Spectrum <- ms_Spectrum[ rank == source, ]
        sp <- get_spectrum_obj(sub_ms_Spectrum, feature_id)
        return(sp)
      })
      sp_diff_sources <- do.call(c, sp_diff_sources)
      return(sp_diff_sources)
    } else{
      stop("method doesn't exist.")
    }
  }else{
    return(NULL)
  }
}

#' Extract MS Spectra of all identified features/peaks
#'
#' @inheritParams plot_Spectra
#'
#' @returns `list` of Spectra objects for all identified features of one sample.
#' @export
exportMSSpectra <- function(features.l, sample_index, type = c("pure", "mixed")[1], method = c("all", "best")[2], max_method = c("contribution", "max_value")[2]) {

  temp_ft <- get_feature_summary(features.l, max_method)
  temp_ft_sample <- temp_ft[ sample == sample_index, ]
  
  MSSpectra <- lapply(unique(temp_ft_sample$featureid), function(featureid) {
    MSSpectrum <- export_featureSpect(features.l, feature_id = featureid, sample_index, type, method, max_method)
  })
  MSSpectra <- do.call(c, MSSpectra)
  return(MSSpectra)
}

# test the three export spectra function
# temp_ft <- get_feature_summary(features.l = features, max_method = "max_value")
# 
# spect <- export_featureSpect(features.l = features, feature_id = "FT1451", sample_index = 1, type = "pure", method = "all", max_method = "contribution")
# ms1_spect <- filterMsLevel(spect, 1L)
# mz(ms1_spect)
# intensity(ms1_spect)
# 
# ms_Spectrum_df <- get_spectra(features.l = features, summary_dt = temp_ft, feature_id = "FT1451",
#                              sample_index = 1,
#                              type = "pure",
#                              method= "best",
#                              max_method = "max_value")
# 
# spect_obj <- get_spectrum_obj(ms_Spectrum_df)
# MsBackendMgf::export(spect, backend = MsBackendMgf(), file = "output.mgf")
# 
# A <- Sys.time()
# all_spectra <- exportMSSpectra(features.l = features, sample_index = 1, type = "pure", method = "best", max_method = "max_value")
# B <- Sys.time()
# B - A
# # Time difference of 7.889102 mins
# some filtering tests on all_spectra
# MsBackendMgf::export(all_spectra, backend = MsBackendMgf(), file = "output.mgf")
# spectraVariables(all_spectra)
# filtered_spectra <- filterValues(all_spectra, "name", "FT2251") # filterValues filter just numerical spectraVariables
# filtered_spect <- all_spectra[spectraData(all_spectra)$name == "FT2251"]
# mz(filtered_spect)
# intensity(filtered_spect)
# ms2_spect <- filterMsLevel(filtered_spect, 2L)
# filterValues(ms2_spect, "isolationWindowTargetMz", 500)

