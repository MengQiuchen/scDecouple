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
#' @param if.plot Logical, indicating whether to plot PC selection line plot or not. Default is FALSE.
#' @param pc.th Numberic, it sets the threshold of x-axis of PC selection line plot. Default is 20.
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
                         top.sdev = 20, top.multimodal = 3, th.sdev = 1, th.multimodal = 0.0125,
                         if.plot=FALSE, pc.th=20) {

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
    if(if.plot){
        colors_ <- c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF','#F39B7FFF','#8491B4FF','#91D1C2FF','#DC0000FF','#7E6148FF','#B09C85FF')
        
        a= Z.pca.sdev[1:pc.th]/sum(Z.pca.sdev)
        b= (Z.pca.mat)[,1:pc.th]%>%apply(2,function(x){x%>%dip.test%>%.[['statistic']]})
        data.df <- cbind(a,b,0)%>%as.data.frame%>%rename_with(~c('Var','Dip','Label'))
        data.df[select.pcs,'Label'] <- 1
        data.df$Dip <- data.df$Dip / 10
        data.df <- data.df%>%rownames_to_column('id')%>%pivot_longer(cols = Var:Dip)%>%mutate(id=gsub('PC','',id))

        data.df$id <- factor(data.df$id, levels=c(1:pc.th))
        data.df$value <- data.df$value * 100
        
        p.pc.percentage <- (ggplot(data.df,aes(x=id,y=value,group=name,color=name))+geom_point(size=0.5)+
            geom_line() +
            scale_color_manual(values = colors_[c(8,4)])+
            scale_y_continuous(
                "Explained Variance (%)", 
                sec.axis = sec_axis(~ . , name = "Modality Score (*0.1)")
              )+
            theme_classic()+
            xlab('PCs'))+
            theme(axis.ticks.y.right = element_line(color =  colors_[8]),
                  axis.text.y.right = element_text(color =  colors_[8]), 
                  axis.title.y.right = element_text(color =  colors_[8]),
                  axis.ticks.y.left = element_line(color = colors_[4]),
                  axis.text.y.left = element_text(color = colors_[4]), 
                  axis.title.y.left = element_text(color =colors_[4])
                 )
        p.pc.percentage%>%print
        
    }


  return(select.pcs)
}
