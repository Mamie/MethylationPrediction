#!/usr/bin/env Rscript
# File name: regressionModel.R
# Author: Mamie Wang
# Date created: 03/01/2018
# Date modified: 03/07/2018

library(glmnet)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
source('../correlation_HM450/computeCorrelation.R')


# Set global options for heatmap plotting
ht_global_opt(heatmap_column_names_gp = gpar(fontsize = 6), 
              heatmap_legend_title_gp= gpar(fontsize = 8),
              heatmap_legend_labels_gp=gpar(fontsize = 6),
              heatmap_row_names_gp = gpar(fontsize = 6))

#' Cluster the methylation matrix with hierachical clustering
#' 
#' @param methylation methylation data matrix 
#' @param distance distance metrics
#' @param method clustering method
#' @param cutoff number of clusters cutoff
#' @return cluster membership of the probes
ClusterMethylation <- function(methylation, distance='euclidean', method='ward.D2', cutoff=20) {
  distances <- dist(methylation, method=distance)
  methylation.hclust <- hclust(distances, method=method)
  cluster <- cutree(methylation.hclust, cutoff)
  return(list(cluster=cluster, hclust=methylation.hclust))
}


#' Average methylation level given clustering
#' 
#' @param methylation methylation data matrix
#' @param cluster a vector indicating cluster membership for the probes
#' @return average methylation level within group
AverageMethylation <- function(methylation, cluster) {
  methylation.grouped <- data.frame(methylation, cluster=cluster) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise_all(funs(mean))
  return(methylation.grouped)
}


#' GLM regression with elastic net regularization of methylation on rnaseq
#' 
#' @param methylation a vector of methylation level
#' @param rnaseq RNAseq data matrix as predictors
#' @param test.idx index of the data that is reserved as testing set
#' @return a glmnet fit object
RegressMethylationOnRNAseq <- function(methylation, rnaseq, test.idx, alpha=1) {
  stopifnot(length(methylation) == dim(rnaseq)[2])
  rnaseq.train <- t(data.matrix(rnaseq[, -test.idx]))
  rnaseq.test <- t(data.matrix(rnaseq[, test.idx]))
  methylation.train <- as.numeric(methylation[-test.idx])
  methylation.test <- as.numeric(methylation[, test.idx])
  glmnet.fit <- cv.glmnet(rnaseq.train, methylation.train, family='gaussian', 
                       alpha=alpha, type.measure='mse')
  return(glmnet.fit)
}


#' Extract non-zero coefficients from a cv.glmnet object
#' 
#' @param glmnet.fit a cv.glmnet object
#' @return a data matrix of nonzero coefficient with row names as gene symbols
ExtractNonzeroCoef <- function(cvglmnet.fit) {
  coef.fit <- data.matrix(predict(cvglmnet.fit, type='coef', s='lambda.min'))
  coef.names <- rownames(coef.fit)[coef.fit != 0]
  coef.nonzero <- coef.fit[coef.fit != 0]
  return(data.frame(gene=coef.names, coef=coef.nonzero))
}


#' Module network visualization of the methylation cluster, averge methylation profile
#' and predictive gene expression level
#' 
#' @param methylation a data matrix of methylation level (each probe is a column)
#' @param avemethyl a column vector of average methylation level across patient
#' @param rnaseq a data matrix of predictive gene RNAseq level (each probe is a column)
#' @param filename the path to which the image file (eps) will be saved
#' @return None (side effect: image saved at filename)
ModuleHeatmap <- function(methylation, avemethyl, rnaseq, filename) {
  rownames(methylation) <- NULL
  rownames(rnaseq) <- NULL
  rownames(avemethyl) <- NULL
  ht1 <- Heatmap(methylation, col=viridis(256), name='methylation', cluster_row=T, 
                 column_title='methylation')
  ha1 <- avemethyl
  ht2 <- Heatmap(ha1, col=viridis(256), width=unit(0.5, 'cm'), name='Average methylation')
  GEP <- rnaseq #scale(rnaseq, center=TRUE, scale=FALSE)
  ht3 <- Heatmap(GEP, col=viridis(256), name='GEP', column_title='GEP', cluster_columns=F)
  setEPS()
  postscript(file=filename)
    draw(ht1 + ht2 + ht3)
  dev.off()
  print(paste('Image saved as', filename))
}


#' Run the regression on all the clusters and plot the module network visulization 
#' 
#' @param methylation a data matrix of methylation level (same number of columns as rnaseq)
#' @param rnaseq a data matrix of RNAseq level
#' @param imagefolder the folder to which the module network image will be saved
#' @param datafolder the folder to which the glmnet object and data will be saved
#' @param distance distance metrics
#' @param method clustering method
#' @param cutoff number of clusters cutoff
#' @param percent.test percent of data to be reserved as testing set
#' @return None (side effect: data and glmnet fit object saved at data folder, 
#' module network images saved at imagefolder)
RunModel <- function(methylation, rnaseq, imagefolder, datafolder, convert2M=F,
                     subsetProbes=NULL, distance='euclidean', method='ward.D2', 
                     cutoff=20, percent.test=0.3, alpha=1, seed=1000) {
  set.seed(seed)
  rnaseq.processed <- PreprocessRNAseq(rnaseq)
  rownames(rnaseq.processed$rnaseq) <- rnaseq.processed$g.geneid[,2]
  methylation.processed <- PreprocessMethylation(methylation, convert2M=convert2M,
                                                 subsetProbes=subsetProbes)
  processed <- matchColumns(rnaseq.processed$rnaseq, 
                            methylation.processed$methylation)
  methylation <- processed$methylation
  rnaseq <- processed$rnaseq
  test.idx <- sample(seq(dim(methylation)[2]), round(dim(methylation)[2]*percent.test))
  clustering <- ClusterMethylation(methylation, distance=distance, method=method, cutoff=cutoff)
  cluster <- clustering$cluster
  avemethyl <- AverageMethylation(methylation, cluster)
  for (i in seq(cutoff)) {
    base.name <- paste0('cluster', i)
    print(paste('Running model for cluster', i))
    path.model <- paste0(datafolder, '/', base.name, '.RData')
    path.fig <- paste0(imagefolder, '/', base.name, '.eps')
    model <- RegressMethylationOnRNAseq(avemethyl[i, -1], rnaseq, test.idx)
    coef <- ExtractNonzeroCoef(model)
    ordering <- order(abs(coef[-1,2]), decreasing=T)
    ordering.genes <- as.character(coef[-1,1][ordering])
    ModuleHeatmap(t(data.matrix(methylation[cluster==i, -test.idx])), 
                  t(data.matrix(avemethyl[i, -1][,-test.idx])), 
                  t(rnaseq[ordering.genes, -test.idx]), 
                  path.fig)
    
    save(model, coef, base.name, file=path.model)
    print(paste('Model saved at', path.model))
  }
  save(test.idx, clustering, avemethyl, 
       file=paste0(datafolder, '/', 'clusterInfo.RData'))
  print(paste('Cluster information saved at', path.model))
}