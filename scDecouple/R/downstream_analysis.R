#' Perform downstream analysis to assess the cellular response of each gene and perform GO enrichment.
#'
#' This function assesses the cellular response of each gene and performs Gene Ontology (GO) enrichment analysis.
#'
#' @param Z.pca.rotation The matrix of variable loadings containing eigenvectors, gene-by-PC, e.g., "Z.pca.results$rotation" from "data_preprocessing" function.
#' @param beta.PC.scDecouple The cellular response of each PC, e.g., "beta.PC.scDecouple" from "scDecouple" function.
#' @param Z.norm.log.sub The log-normalized control group gene expression matrix, cell-by-gene, e.g., "Z.norm.log.sub" from "data_preprocessing" function.
#' @param Y.norm.log.sub The log-normalized perturbation group gene expression matrix, cell-by-gene, e.g., "Y.norm.log.sub" from "data_preprocessing" function.
#' @param genes.variable A list of variable genes selected for analysis.
#' @param degenes.num Optional. The number of genes to be included in GO enrichment analysis. Default is 1500.
#'
#' @return A list containing the following results:
#'   \itemize{
#'     \item gene.rank: A vector of genes sorted by their response to perturbation.
#'     \item GO.enrich: GO enrichment results using the enrichr function in the enrichR library.
#'   }
#'
#' @import tidyverse
#' @import enrichR
#'
#' @examples
#' # Example usage:
#' results <- downstream_analysis(Z.pca.rotation, beta.PC.scDecouple, Z.norm.log.sub, Y.norm.log.sub, genes.variable, degenes.num)
#'
#' @export
downstream_analysis <- function(Z.pca.rotation, beta.PC.scDecouple, Z.norm.log.sub, Y.norm.log.sub, genes.variable,
                                degenes.num = 1500) {
  
  library(tidyverse)
  
  # Calculate the response of each gene by back projection to gene space
  beta.variable.scDecouple <- beta.PC.scDecouple %*% t(Z.pca.rotation) - colMeans(Z.norm.log.sub)
  beta.gene.fc <- colMeans(Y.norm.log.sub) - colMeans(Z.norm.log.sub)
  
  beta.gene.scDecouple <- beta.gene.fc
  beta.gene.scDecouple[genes.variable] <- beta.variable.scDecouple[, genes.variable]
  
  # Gene ranking
  results.genesort <- beta.gene.scDecouple[beta.gene.scDecouple %>% abs %>% sort(decreasing = TRUE) %>% names]
  results.genesort.fc <- beta.gene.fc[beta.gene.fc %>% abs %>% sort(decreasing = TRUE) %>% names]
  
  # Gene enrichment using enrichR
  library('enrichR')
  suppressMessages({
    GO.enrich <- results.genesort[1:degenes.num] %>% names %>% enrichr('GO_Biological_Process_2018') %>% .$GO_Biological_Process_2018
  })
  
  return(list(gene.rank = results.genesort.fc,
              GO.enrich = GO.enrich))
}
