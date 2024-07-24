# Rank of factorization methods

#' determine the rank of factorization using concordance metric
#'
#' @param mat.m matrix
#'
#' @import reticulate
#' @return rank of mat.m
#'
concordance_rank <- function(mat.m){
  
  mat.m <- t(apply(mat.m, 1, function(x) x / max(x)));
  np_data_matrix <- r_to_py(mat.m);
  my_comp_max <- min( as.integer(ncol(mat.m)), as.integer(nrow(mat.m)) );
  res <- perform_rank_determination(data = np_data_matrix, comp_min = as.integer(2), comp_max = my_comp_max,
                                    iter_max = as.integer(20), n_runs = as.integer(10), method = "concordance", save_res = FALSE);
  rank <- res[[4]];
  # if(rank == 999) { rank = 2 };
  
  return(rank)
}

