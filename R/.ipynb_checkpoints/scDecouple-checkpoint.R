#' Perform scDecouple analysis to deconvolute cellular response from infection bias.
#'
#' This function uses scDecouple to deconvolute the cellular response from infection bias in two cell-by-gene matrices.
#'
#' @param Z.pca.mat The PC matrix of the control group, cell-by-PC, e.g., result from "data_preprocessing" function.
#' @param Z.pca.rotation The matrix of variable loadings containing eigenvectors, gene-by-PC, e.g., "Z.pca.results$rotation" from "data_preprocessing" function.
#' @param Y.norm.log.sub The PC matrix of the perturbation group, cell-by-PC.
#' @param select.pcs A character vector containing the names of selected PCs.
#' @param seed_ Optional. The random seed for the EM algorithm. Default is 0.
#' @param if.plot Logical, indicating whether to plot scatter plot or not. Default is FALSE.
#'
#' @return A list containing the following results:
#'   \itemize{
#'     \item ratio.changes: The ratio changes of the perturbation group compared to the control group for each cluster.
#'     \item beta.PC.scDecouple: The cellular response of each PC to the perturbation.
#'   }
#'
#' @import tidyverse
#' @import mixtools
#'
#' @examples
#' # Example usage:
#' result <- scDecouple(Z.pca.mat, Z.pca.rotation, Y.norm.log.sub, select.pcs, seed_)
#'
#' @export
scDecouple <- function(Z.pca.mat, Z.pca.rotation, Y.norm.log.sub, select.pcs, seed_ = 0, if.plot=FALSE) {
  
  library(tidyverse)
  
  # Projection
  Z.PCA <- Z.pca.mat
  
  # Map Y to PC space
  Y.PCA <- Y.norm.log.sub %*% Z.pca.rotation
  
  library(mixtools)
  
  # Control group cluster properties estimation
  set.seed(seed_)
  Z.EM <- mvnormalmixEM(Z.PCA[, select.pcs], k = 2)
  
  z.labels <- Z.EM$posterior %>% apply(1, which.max)
  
  # Deconvolution
  beta.PC <- colMeans(Y.PCA) - colMeans(Z.PCA)
  
  Y.EM <- get_beta_EM_sep_method(y = Y.PCA[, select.pcs], number.component = 2,
                                 em = Z.EM$mu, es = Z.EM$sigma, el = Z.EM$lambda,
                                 beta_init = beta.PC[select.pcs])
  beta.scDecouple <- Y.EM$beta
  
  # Show results
  # Ratio changes
  ratio.changes <- Y.EM$el - Z.EM$lambda
  
  # Cellular responses on PCs
  beta.PC.scDecouple <- beta.PC
  beta.PC.scDecouple[select.pcs] <- beta.scDecouple
  
  # plots of control cells
  if(if.plot){
        
        mat.df.2d <- cbind(Z.PCA[,select.pcs[c(1,2)]],(z.labels))%>%
                as.data.frame%>%rename_with(~c('PC_selected_1','PC_selected_2','label'))
        p.pcs <- ggplot(mat.df.2d,aes(x=PC_selected_1,y=PC_selected_2,color=as.character(z.labels)))+
                geom_point(cex=0.1)+   
                scale_color_manual(values=colors_[c(4,8)])+
                scale_fill_manual(values=colors_[c(4,8)])+
                theme_classic()+
                #xlab('PC1 value')+ylab('Frequency')+
                theme(aspect.ratio = 1,legend.position = 'none',
                     axis.text=element_text(size=8,color = 'black'),
                     axis.title=element_text(size=8,face='bold',color = 'black'))
        p.pcs%>%print
        
    }
    
  return(list(ratio.changes = ratio.changes,
              beta.PC.scDecouple = beta.PC.scDecouple))
}
