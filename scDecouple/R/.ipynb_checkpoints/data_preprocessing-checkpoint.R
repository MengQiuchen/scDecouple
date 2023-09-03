#' Data preprocessing for scDecouple analysis.
#'
#' This function performs data preprocessing for the scDecouple analysis on two cell-by-gene matrices, Y (perturbation group) and Z (control group).
#'
#' @param Y The perturbation group cell-by-gene matrix.
#' @param Z The control group cell-by-gene matrix.
#'
#' @return A list containing the following preprocessed data:
#'   \itemize{
#'     \item Z.norm.log.sub: Normalized and log-transformed data matrix for the control group with variable genes.
#'     \item Y.norm.log.sub: Normalized and log-transformed data matrix for the perturbation group with variable genes.
#'     \item genes.variable: List of variable genes selected for analysis.
#'     \item Z.pca.results: PCA results for the control group data.
#'   }
#'
#' @examples
#' # Example usage:
#' preprocessed_data <- data_preprocessing(Y, Z)
#'
#' @import Seurat
#' @export
data_preprocessing <- function(Y, Z) {

  # Step 1: Select Variable Genes using Seurat
  library(Seurat)
  suppressWarnings({
    Z.seurat <- CreateSeuratObject(counts = t(Z), min.cells = 0, min.features = 0)
    Z.seurat <- NormalizeData(Z.seurat)
    Z.seurat <- FindVariableFeatures(Z.seurat, selection.method = "vst", nfeatures = 600)
    genes.variable <- VariableFeatures(Z.seurat)
  })

  # Step 2: Normalize and Log-Transform Data
  Z.norm <- (Z / rowSums(Z) * 10000)
  Y.norm <- (Y / rowSums(Y) * 10000)
  Z.norm.log <- log(Z.norm + 1)
  Y.norm.log <- log(Y.norm + 1)

  # Step 3: Subset Data with Variable Genes
  Z.norm.log.sub <- Z.norm.log[, genes.variable]
  Y.norm.log.sub <- Y.norm.log[, genes.variable]

  # Step 4: Perform PCA on Control Group Data
  Z.pca.results <- prcomp(Z.norm.log.sub, center = TRUE, scale. = FALSE)

  return(list(Z.norm.log.sub = Z.norm.log.sub,
              Y.norm.log.sub = Y.norm.log.sub,
              genes.variable = genes.variable,
              Z.pca.results = Z.pca.results))
}
