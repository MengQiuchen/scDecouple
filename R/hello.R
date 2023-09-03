suppressPackageStartupMessages({
  library(tidyverse)
  library(enrichR)
  library(mixtools)
  library(diptest)
  library(Seurat)
})

# decoupling process
get_beta_EM_sep_method <- function(y, number.component, em, es, el, beta_init,
                                   max.step = 100, fix_es = FALSE){

  library(mixtools)
    
  if(!length(dim(y))){# one-dim
    number.sample <- length(y)
    dim <- 1
  }else{
    number.sample <- dim(y)[1]
    dim <- dim(y)[2]
  }

  Q_before <- 0
  beta <- beta_init
  ey <- matrix(NA, nrow = number.sample, ncol = number.component)
  dy <- matrix(NA, nrow = number.sample, ncol = number.component)
  for (k in 1:max.step) {

    for (j in 1:number.component) {
      if(dim==1){
        dy[, j] <- (dnorm(data.matrix(y), as.matrix(em[[j]] +
                                                      beta), es[[j]])) * el[j]
      }else{
        dy[, j] <- (dmvnorm(data.matrix(y), as.matrix(em[[j]] +
                                                        beta), es[[j]])) * el[j]
      }

    }
    for (i in 1:number.sample) {
      sy <- sum(dy[i, ])
      if (sy == 0) {
        sy <- 1e-308
      }
      ey[i, ] <- dy[i, ]/sy
    }
    em_sum <- 0


    flag=TRUE
    for (j in 1:number.component) {
      el[j] <- sum(ey[, j])/number.sample

      if(!sum(ey[,j])){
        flag=FALSE
        print('WARNING:NOT CONVERGENT!')
        break
      }
      if(dim==1){
        em_sum <- em_sum + sum(ey[, j] * y)/sum(ey[,
                                                   j])
      }else{
        em_sum <- em_sum + colSums(ey[, j] * y)/sum(ey[,
                                                       j])
      }

    }
    if(!flag){ break }
    beta <- em_sum/number.component - colMeans(matrix(unlist(em),
                                                      ncol = dim, byrow = TRUE))
    if (!fix_es) {
      for (j in 1:number.component) {
        es_sum <- matrix(0, nrow = dim, ncol = dim)
        for (i in 1:number.sample) {
          if(dim==1){
            es_sum <- es_sum + ey[i, j] * (t(y[i] - t(as.matrix(em[[j]] +
                                                                  beta))) %*% (y[i] - t(as.matrix(em[[j]] +
                                                                                                    beta))))
          }else{
            es_sum <- es_sum + ey[i, j] * (t(y[i, ] - t(as.matrix(em[[j]] +
                                                                    beta))) %*% (y[i, ] - t(as.matrix(em[[j]] +
                                                                                                        beta))))
          }

        }
        es[[j]] <- es_sum/sum(ey[, j])
      }
    }
    Q_after <- sum(log(rowSums(dy)))
    if (abs(Q_after - Q_before) < 1e-06) {
      break
    }
    else {
      Q_before <- Q_after
    }
  }
  result <- list(beta = beta, es = es, el = el, em = em, Q_after = Q_after,
                 ey = ey)
}

# step1: data preprocessing

# the input is
## Y: the perturbation group, cell-by-gene matrix
## Z: the control group, cell-by-gene matrix
data_preprocessing <- function(Y,Z){

  # the code shows the steps of preprocessing
  ## 1. normalization and transformation
  ## 2. variable gene selection

  ## the following code shows a demo of data preprocessing

  # select variable genes
  ## using Seurat to do preprocessing
    
    library(Seurat)
  suppressWarnings({

    Z.seurat <- CreateSeuratObject(counts = t(Z), min.cells = 0, min.features = 0)
    Z.seurat <- NormalizeData(Z.seurat)
    Z.seurat <- FindVariableFeatures(Z.seurat, selection.method = "vst", nfeatures = 600)

    genes.variable <- VariableFeatures(Z.seurat)
  })

  # ## one method to do normalization & transfromation
  # Z.norm.log <- as.matrix(Z.seurat@assays$RNA@data)

  # ## another method
  #### normalization & transformation
  ## normalization another method to do normalization
  Z.norm <- (Z/rowSums(Z)*10000)
  Y.norm <- (Y/rowSums(Y)*10000)
  ## log transformation
  Z.norm.log <- log(Z.norm+1)
  Y.norm.log <- log(Y.norm+1)

  ## subset data matrix with variable genes
  Z.norm.log.sub <- Z.norm.log[,genes.variable]
  Y.norm.log.sub <- Y.norm.log[,genes.variable]

  ## PCA
  Z.pca.results <- prcomp(Z.norm.log.sub, center=TRUE, scale. = FALSE)

  return(list(Z.norm.log.sub=Z.norm.log.sub,
              Y.norm.log.sub=Y.norm.log.sub,
              genes.variable=genes.variable,
              Z.pca.results=Z.pca.results))

}


# step2: PC select

## the input is:
## 1. Z.pca.mat: PC matrix of control group, cells-by-PCs
##### for example, the "Z.pca.results$x" from "data_preprocessing" function
## 2.(optional) Z.pca.sdev: the standard deviation of PCs, the default is NULL
##### for example, the "Z.pca.results$sdev" from "data_preprocessing" function
##### if is not provided, we will calculated from Z.pca.mat
## 3.(optional) using_rank: using rank (select PCs with top sdev &multimodality) or using value (set the threshold for sdev &mutlimodality) to select PCs, the default is TRUE
## 4.(optinal) top.sdev: when using_rank=TRUE, to set the threshold of standard deviations of PCs, the default is 20
## 5.(optinal) top.multimodal: when using_rank=TRUE, to set the threshold of dip statistics, the default is 3
## 6.(optional) th.sdev: when using_rank=FALSE, to set the threshold of standard deviations of PCs, the default is 1
## 7.(optional) th.multimodal: when using_rank=FALSE, to set the threshold of dip statistics, the default is 0.0125

## the code shows the steps of PC selection based on explained variance and multimodality

pc_selection <- function(Z.pca.mat,
                         Z.pca.sdev=NULL,using_rank=TRUE,
                         top.sdev=20, top.multimodal=3, th.sdev=1, th.multimodal=0.0125){

    library(diptest)
    library(tidyverse)
  # if Z.pca.sdev is not provided, we calculated by Z.pca.mat
  if(!length(Z.pca.sdev)){
    Z.pca.sdev <- Z.pca.mat%>%apply(2,sd)%>%as.numeric
  }

  if(using_rank){
    # use rank to select PCs
    library(diptest)
    select.pcs <- ((Z.pca.mat)[,1:top.sdev]%>%apply(2,function(x){x%>%dip.test%>%.[['statistic']]})%>%
                     sort(decreasing = TRUE))[1:top.multimodal]%>%names

  }else{
    # use value to select PCs
    # selecte PCs by sdev
    pc.select.sdev <- which(Z.pca.sdev > th.sdev)
    # selecte PCs by multimodality
    library(diptest)
    select.pcs <- which((Z.pca.mat)[,pc.select.sdev]%>%
                          apply(2,function(x){x%>%dip.test%>%.[['statistic']]})> th.multimodal)%>%names
  }

  return(select.pcs)


}


# step3: deconvolution

# this part used scDecouple to deconvolute the cellular response from infection bias

## input:
## 1. Z.pca.mat: PC matrix of control group, cell-by-PC.
##### for example, the "Z.pca.results$x" from "data_preprocessing" function
## 2. Z.pca.rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors), gene-by-PC.
#### for example, the "Z.pca.results$rotation" from "data_preprocessing" function, or "rotation" element from "prcomp" function.
## 3. Y.norm.log.sub: PC matrix of perturbation group, cell-by-PC.
## 4. seed_:(optional) the random seed for EM algorithm. Default is 0.

## output:
## 1. ratio.changes: the ratio changes of perturbation group compared with control group on chluster 1 to cluster n.
## 2. beta.PC.scDecouple: the cellular response of each PC to the perturbation. the results were stored in list of double.

scDecouple <- function(Z.pca.mat, Z.pca.rotation,Y.norm.log.sub,select.pcs, seed_=0){
    library(tidyverse)
  ## projection
  Z.PCA <- Z.pca.mat

  ## map Y to PC space
  Y.PCA <- Y.norm.log.sub  %*% Z.pca.rotation

  library(mixtools)

  ## control group cluster proporties estimation
  set.seed(seed_)
  Z.EM <- mvnormalmixEM(Z.PCA[,select.pcs],k = 2)

  z.labels <- Z.EM$posterior%>%apply(1,which.max)

  # plots of control cells

  # deconvolutio
  beta.PC <- colMeans(Y.PCA) - colMeans(Z.PCA)

  Y.EM <- get_beta_EM_sep_method(y = Y.PCA[,select.pcs],number.component = 2,
                                 em = Z.EM$mu, es = Z.EM$sigma, el = Z.EM$lambda,
                                 beta_init = beta.PC[select.pcs] )
  beta.scDecouple <- Y.EM$beta

  # show results
  ## ratio.changes
  ratio.changes <- Y.EM$el - Z.EM$lambda

  ## cellular responses on PCs
  beta.PC.scDecouple <- beta.PC
  beta.PC.scDecouple[select.pcs] <- beta.scDecouple

  return(list(ratio.changes=ratio.changes,
              beta.PC.scDecouple=beta.PC.scDecouple
  ))
}

# step4: downstream analysis

## cellular response of each genes
## input:
## 1. Z.pca.rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors), gene-by-PC.
#### for example, the "Z.pca.results$rotation" from "data_preprocessing" function, or "rotation" element from "prcomp" function.
## 2. beta.PC.scDecouple: the cellular response of each PC
#### for example, the "beta.PC.scDecouple" from "scDecouple" function
## 3. Z.norm.log.sub: the log-normalized control group gene expression matrix, cell-by-gene.
#### for example, the "Z.norm.log.sub" from "data_preprocessing" function
## 4. Y.norm.log.sub: the log-normalized perturbation group gene expression matrix, cell-by-gene.
#### for example, the "Y.norm.log.sub" from "data_preprocessing" function
## 5. degenes.num: (optional) the number of degenes to calculated GO enrichment. (default is 1500)

## output:
## 1. gene.rank: the rank of genes sorted by their response to perturbation
## 2. GO.enrich: the enrichment results provided by enrichr function in enrichR library.
downstream_analysis <- function(Z.pca.rotation, beta.PC.scDecouple, Z.norm.log.sub, Y.norm.log.sub,genes.variable,
                                degenes.num=1500){

    library(tidyverse)
  # calculate the response of each gene by back projection to gene space
  beta.variable.scDecouple <- beta.PC.scDecouple %*% t(Z.pca.rotation) - colMeans(Z.norm.log.sub)
  beta.gene.fc <- colMeans(Y.norm.log.sub) - colMeans(Z.norm.log.sub)

  beta.gene.scDecouple <- beta.gene.fc
  beta.gene.scDecouple[genes.variable] <- beta.variable.scDecouple[,genes.variable]

  ## gene ranking
  results.genesort <- beta.gene.scDecouple[beta.gene.scDecouple%>%abs %>% sort(decreasing = TRUE)%>%names]

  results.genesort.fc <- beta.gene.fc[beta.gene.fc%>%abs %>% sort(decreasing = TRUE)%>%names]

  ## ## gene enrichment
  library('enrichR')
  degenes.num <- 1500
  suppressMessages({
    GO.enrich <- results.genesort[1:degenes.num]%>%names%>%enrichr('GO_Biological_Process_2018')%>%.$GO_Biological_Process_2018
  })

  return(list(gene.rank = results.genesort.fc,
              GO.enrich = GO.enrich))
}




## Y: the perturbation group, cell-by-gene matrix
## Z: the control group, cell-by-gene matrix
## 3.(optional) using_rank: using rank (select PCs with top sdev &multimodality) or using value (set the threshold for sdev &mutlimodality) to select PCs, the default is TRUE
## 4.(optinal) top.sdev: when using_rank=TRUE, to set the threshold of standard deviations of PCs, the default is 20
## 5.(optinal) top.multimodal: when using_rank=TRUE, to set the threshold of dip statistics, the default is 3
## 6.(optional) th.sdev: when using_rank=FALSE, to set the threshold of standard deviations of PCs, the default is 1
## 7.(optional) th.multimodal: when using_rank=FALSE, to set the threshold of dip statistics, the default is 0.0125
## 5. degenes.num: (optional) the number of degenes to calculated GO enrichment. (default is 1500)
## 4. seed_:(optional) the random seed for EM algorithm. Default is 0.

run_scDecouple <- function(Y,Z,
                           using_rank=TRUE,top.sdev=30,top.multimodal=3,th.sdev=1,
                           th.multimodal=0.0125,degenes.num=1500,seed_=0){

  res.step1 <- data_preprocessing(Y,Z)
  res.step2 <- pc_selection(Z.pca.mat = res.step1$Z.pca.results$x,Z.pca.sdev = res.step1$Z.pca.results$sdev,using_rank = TRUE)
  res.step3 <- scDecouple(Z.pca.mat = res.step1$Z.pca.results$x, Z.pca.rotation = res.step1$Z.pca.results$rotation,Y.norm.log.sub = res.step1$Y.norm.log.sub,select.pcs=res.step2)
  res.step4 <- downstream_analysis(Z.pca.rotation =res.step1$Z.pca.results$rotation,beta.PC.scDecouple = res.step3$beta.PC.scDecouple,Z.norm.log.sub = res.step1$Z.norm.log.sub,Y.norm.log.sub =res.step1$Y.norm.log.sub,genes.variable=res.step1$genes.variable)

  return(list(processed_data=res.step1,
              selected_pcs=res.step2,
              decouple_results=res.step3,
              downstream_analysis=res.step4))
}
