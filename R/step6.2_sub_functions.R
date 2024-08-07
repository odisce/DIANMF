#' Delete rows that have only one non-zero element
#'
#' @param mat matrix
#'
#' @return filtered matrix
filter_rows <- function(mat) {
  rows_to_keep <- apply(mat, 1, function(row) sum(row != 0)) > 2
  mat[rows_to_keep, , drop = FALSE]  # Use drop = FALSE to maintain matrix structure
}


#' Delete rows that have at least 3 consecutive non-zero values
#'
#' @param row 
#'
#' @return filtered row
#' @export
has_four_consecutive_non_zero <- function(row) {
  # Find consecutive non-zero values
  non_zero_streaks <- rle(row != 0)
  any(non_zero_streaks$lengths[non_zero_streaks$values] >= 4)
}
# mat <- matrix(c(
#   0, 1, 2, 3, 0,
#   2, 0, 0, 0, 0,
#   0, 3, 4, 0, 0,
#   1, 0, 2, 0, 3
# ), nrow = 4, byrow = TRUE)
# row_filter <- apply(mat, 1, has_four_consecutive_non_zero)
# filtered_mat <- mat[row_filter, , drop = FALSE]


#' filter matrices in a list
#'
#' @param mat 
#'
#' @return filtered matrix
filter_matrix_rows <- function(mat) {
  row_filter <- apply(mat, 1, has_four_consecutive_non_zero)
  mat[row_filter, , drop = FALSE]
}


#' calculate the elution profiles correlation to choose the good component
#'
#' @param chromo_main 
#' @param chromos 
#'
#' @return index of the corresponding compound
#' @importFrom stats cor
elutions_corelation <- function(chromo_main, chromos){
  correlations <- c()
  for (i in 1:nrow(chromos)) {
    correlations[i] <- stats::cor(chromo_main, chromos[i, ])
  }
  comp <- which.max(correlations)
  return(comp)
}

