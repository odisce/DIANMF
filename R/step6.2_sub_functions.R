#' calculate the elution profiles correlation to choose the good component
#'
#' @param chromo_main numeric_vector MS1 good elution profile
#' @param chromos numeric_vectorMS2 elution profiles
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

