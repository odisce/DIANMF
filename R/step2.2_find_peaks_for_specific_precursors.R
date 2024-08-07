#' load the compounds file as data-frames/table
#' 
#' @param compound.c character path to load the spiked compounds table
#' 
#' @return data-frame/table spiked compounds table
#' @export
#' @importFrom utils read.table
load_compounds_function.tb <- function(compound.c){
  m <- read.table(compound.c,
                  sep = "\t",
                  quote = "",
                  header = TRUE,
                  row.names = 1)
  m[, "Retention.time.s."] <- m[, "Retention.time..min."] * 60
  rownames(m) <- make.names(rownames(m), unique = TRUE)
  
  return(m)
}

#-------------------------------------------------------------------------------

#' create list for the spiked precursors/compounds
#'
#' @param compound.df table, data-frame of the spiked compounds
#' 
#' @return empty list of all precursors I'm searching for
#' @export
create_list_function.l <- function(compound.df){
  compound_names.vc <- paste0(rownames(compound.df),
                              "_",
                              compound.df[, "Composition"],
                              "_",
                              compound.df[, "Precursor.mz"],
                              "_",
                              compound.df[, "Retention.time.s."])
  compound.l <- vector(mode = "list", length = length(compound_names.vc))
  names(compound.l) <- compound_names.vc
  
  return(compound.l)
}

#-------------------------------------------------------------------------------

#' search for the peaks suggestions of every compound using the mz and retention time info
#'
#' @param compounds.l list empty, every level related to a compounds
#' @param matrix.n matrix of peaks
#' @param mz_tol.n numeric mz tolerance
#' @param  rt_tol.n numeric rt tolerance
#' 
#' @return filled list of peak suggestions for all levels
#' @export
search_for_ROIs_function.l <- function(compounds.l, matrix.n, rt_tol.n){   # I want to make this function more critical, when we have several match
  
  for (compound.c in names(compounds.l)) {
    precursor_info.vc <- unlist(strsplit(compound.c, split = "_"))
    mz.n <- as.numeric(precursor_info.vc[3])
    rt.n <- as.numeric(precursor_info.vc[4])
    
    # mz_tol.n <- mz.n * 1e-5
    mz_tol.n <- 0.05
      
    compounds.l[[compound.c]] <- which(matrix.n[, "mzmin"] - mz_tol.n <= mz.n &
                                         matrix.n[, "mzmax"]+ mz_tol.n >=  mz.n &
                                         matrix.n[, "rtmin"] - rt_tol.n <= rt.n  &
                                         matrix.n[, "rtmax"] + rt_tol.n >= rt.n)
  }
  return(compounds.l)
}

#-------------------------------------------------------------------------------