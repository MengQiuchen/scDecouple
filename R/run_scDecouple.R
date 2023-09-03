#' Run the scDecouple analysis on two cell-by-gene matrices.
#'
#' This function performs the scDecouple analysis on two cell-by-gene matrices, Y (perturbation group) and Z (control group).
#'
#' @param Y The perturbation group cell-by-gene matrix.
#' @param Z The control group cell-by-gene matrix.
#' @param using_rank Logical, indicating whether to use rank-based PC selection or not. Default is TRUE.
#' @param top.sdev Integer, the threshold of standard deviations of PCs when using_rank is TRUE. Default is 30.
#' @param top.multimodal Numeric, the threshold of dip statistics when using_rank is TRUE. Default is 3.
#' @param th.sdev Numeric, the threshold of standard deviations of PCs when using_rank is FALSE. Default is 1.
#' @param th.multimodal Numeric, the threshold of dip statistics when using_rank is FALSE. Default is 0.0125.
#' @param degenes.num Numeric, the number of degenes to calculate GO enrichment. Default is 1500.
#' @param seed_ Numeric, the random seed for the EM algorithm. Default is 0.
#'
#' @return A list containing the following elements:
#'   \itemize{
#'     \item processed_data: Results of data preprocessing.
#'     \item selected_pcs: Selected principal components.
#'     \item decouple_results: Results of the scDecouple analysis.
#'     \item downstream_analysis: Results of downstream analysis.
#'   }
#'
#' @examples
#' # Example usage:
#' result <- run_scDecouple(Y, Z)
#'
#' @export
run_scDecouple <- function(Y, Z,
                           using_rank = TRUE, top.sdev = 30, top.multimodal = 3, th.sdev = 1,
                           th.multimodal = 0.0125, degenes.num = 1500, seed_ = 0) {
  
  # Step 1: Data Preprocessing
  res.step1 <- data_preprocessing(Y, Z)
  
  # Step 2: PC Selection
  res.step2 <- pc_selection(Z.pca.mat = res.step1$Z.pca.results$x,
                            Z.pca.sdev = res.step1$Z.pca.results$sdev,
                            using_rank = using_rank,
                            top.sdev = top.sdev,
                            top.multimodal = top.multimodal)
  
  # Step 3: scDecouple Analysis
  res.step3 <- scDecouple(Z.pca.mat = res.step1$Z.pca.results$x,
                          Z.pca.rotation = res.step1$Z.pca.results$rotation,
                          Y.norm.log.sub = res.step1$Y.norm.log.sub,
                          select.pcs = res.step2)
  
  # Step 4: Downstream Analysis
  res.step4 <- downstream_analysis(Z.pca.rotation = res.step1$Z.pca.results$rotation,
                                   beta.PC.scDecouple = res.step3$beta.PC.scDecouple,
                                   Z.norm.log.sub = res.step1$Z.norm.log.sub,
                                   Y.norm.log.sub = res.step1$Y.norm.log.sub,
                                   genes.variable = res.step1$genes.variable)
  
  return(list(processed_data = res.step1,
              selected_pcs = res.step2,
              decouple_results = res.step3,
              downstream_analysis = res.step4))
}
