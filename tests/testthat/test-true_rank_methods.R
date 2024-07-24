test_that("rank method test", {
  
  # install.packages("pracma")
  # library(pracma)
  # 
  # m <- matrix(c(1,0,0,0,1,1,1,1,1), nrow = 3)
  # pracma::Rank(m)
  # #
  # 
  # Matrix::rankMatrix(m, method = 'maybeGrad')
  # 
  
  
  # mnist_dataSet <- read_mnist(
  #   path = NULL,
  #   download = FALSE,
  #   destdir = tempdir(),
  #   url = "https://www2.harvardx.harvard.edu/courses/IDS_08_v2_03/",
  #   keep.files = TRUE
  # )
  # # Note that the actual rank for this data set is 10.
  # 
  # mat.m <- mnist_dataSet$test$images
  # mat.m <- t(apply(mat.m, 1, function(x) x / max(x)));
  # np_data_matrix <- r_to_py(mat.m);
  # my_comp_max <- min( as.integer(ncol(mat.m)), as.integer(nrow(mat.m)) );
  # res <- perform_rank_determination(data = np_data_matrix, comp_min = as.integer(2), comp_max = my_comp_max,
  #                                   iter_max = as.integer(20), n_runs = as.integer(10), method = "concordance", save_res = FALSE);
  # rank <- res[[4]];
})










