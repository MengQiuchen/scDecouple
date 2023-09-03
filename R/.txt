#' Perform principal component (PC) selection for scDecouple analysis.
#'
#' This function selects principal components (PCs) based on explained variance and multimodality measures.
#'
#' @param Z.pca.mat The PC matrix of the control group, cells-by-PCs, e.g., result from "data_preprocessing" function.
#' @param Z.pca.sdev Optional. The standard deviation of PCs. If not provided, it will be calculated from Z.pca.mat.
#' @param using_rank Logical, indicating whether to use rank-based PC selection or not. Default is TRUE.
#' @param top.sdev Numeric. When using_rank is TRUE, it sets the threshold of standard deviations of PCs. Default is 20.
#' @param top.multimodal Numeric. When using_rank is TRUE, it sets the threshold of dip statistics. Default is 3.
#' @param th.sdev Numeric. When using_rank is FALSE, it sets the threshold of standard deviations of PCs. Default is 1.
#' @param th.multimodal Numeric. When using_rank is FALSE, it sets the threshold of dip statistics. Default is 0.0125.
#'
#' @return A character vector containing the names of selected PCs.
#'
#' @import diptest
#' @import tidyverse
#'
#' @examples
#' # Example usage:
#' selected_pcs <- pc_selection(Z.pca.mat, Z.pca.sdev, using_rank, top.sdev, top.multimodal, th.sdev, th.multimodal)
#'
#' @export
pc_selection <- function(Z.pca.mat,
                         Z.pca.sdev = NULL, using_rank = TRUE,
                         top.sdev = 20, top.multimodal = 3, th.sdev = 1, th.multimodal = 0.0125) {

  library(diptest)
  library(tidyverse)

  # Calculate Z.pca.sdev if not provided
  if (!length(Z.pca.sdev)) {
    Z.pca.sdev <- Z.pca.mat %>% apply(2, sd) %>% as.numeric
  }

  if (using_rank) {
    # Use rank to select PCs
    select.pcs <- ((Z.pca.mat)[, 1:top.sdev] %>% apply(2, function(x) {
      x %>% dip.test %>% .[['statistic']]
    }) %>%
      sort(decreasing = TRUE))[1:top.multimodal] %>% names

  } else {
    # Use value to select PCs
    # Select PCs by sdev
    pc.select.sdev <- which(Z.pca.sdev > th.sdev)
    # Select PCs by multimodality
    select.pcs <- which((Z.pca.mat)[, pc.select.sdev] %>%
                          apply(2, function(x) {
                            x %>% dip.test %>% .[['statistic']]
                          }) > th.multimodal) %>% names
  }

  return(select.pcs)
}
