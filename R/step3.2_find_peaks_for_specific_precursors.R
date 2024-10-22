#' Load the compounds table.
#' 
#' @param compound.c file path of the spiked compounds table
#' 
#' @return data.frame/data.table spiked compounds table.
#' 
#' @export
#' 
#' @examples
#' compounds_tb <- system.file("extdata", "Barbier_supplementary_tableS1.tsv", package = "DIANMF")
#' 
#' load_compounds_function.tb(compounds_tb)
#' 
#' @importFrom utils read.table
load_compounds_function.tb <- function(compound.c){
  m <- utils::read.table(compound.c,
                  sep = "\t",
                  quote = "",
                  header = TRUE,
                  row.names = 1)
  m[, "Retention.time.s."] <- m[, "Retention.time..min."] * 60
 
  rownames(m) <- make.names(rownames(m), unique = TRUE)
  colnames(m) <- c("composition", "mz", "rt.min", "NCE", "fragment1.mz", "fragment2.mz", "fragment3.mz", "rt.sec")
  
  return(m)
}


#' Search for the peaks suggestions for compounds using the mz and retention time info.
#'
#' @param compounds.df `data.frame` 
#' @inheritParams prepare_ms1_peaks
#' @param mz_tol.n `numeric(1)` mz tolerance.
#' @param rt_tol.n `numeric(1)` retention time tolerance.
#'
#' @return MS1 selected peaks `matrix`.
#' 
#' @export
search_for_ROIs_function.l <- function(compounds.df, ms1_peaks, mz_tol.n = 5, rt_tol.n = 5) {
  res <- lapply(1:nrow(compounds.df), function(i) {
    
    mz.n <- as.numeric(compounds.df[i, "mz"])
    rt.n <- as.numeric(compounds.df[i, "rt.sec"])
    
    matches <- ms1_peaks[
      (ms1_peaks[, "mzmin"] - mz_tol.n <= mz.n & ms1_peaks[, "mzmax"] + mz_tol.n >= mz.n) &
        (ms1_peaks[, "rtmin"] - rt_tol.n <= rt.n & ms1_peaks[, "rtmax"] + rt_tol.n >= rt.n), 
      , drop = FALSE]
    
    if (nrow(matches) == 0) {
      return(NULL)  # No matches found
    } else if (nrow(matches) == 1) {
      rownames(matches) <- NULL
      return(matches)  # Only one match, return it
    } else {
      mz_centers <- rowMeans(matches[, c("mzmin", "mzmax")])
      closest_index <- which.min(abs(mz_centers - mz.n))
      
      selected_match <- matches[closest_index, , drop = FALSE]
      rownames(selected_match) <- NULL
      return(selected_match)
    }
  })
  
  res <- do.call(rbind, res)
  
  return(res)
}
