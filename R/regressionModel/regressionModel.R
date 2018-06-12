#!/usr/bin/env Rscript
# File name: regressionModel.R
# Author: Mamie Wang
# Date created: 03/01/2018
# Date modified: 03/07/2018

library(glmnet)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(vbsr)
library(qvalue)
#source('../correlation_HM450/computeCorrelation.R')


# Set global options for heatmap plotting
ht_global_opt(heatmap_column_names_gp = gpar(fontsize = 2), 
              heatmap_legend_title_gp= gpar(fontsize = 8),
              heatmap_legend_labels_gp=gpar(fontsize = 6),
              heatmap_row_names_gp = gpar(fontsize = 2))

#' Cluster the methylation matrix with hierachical clustering
#' 
#' @param methylation methylation data matrix 
#' @param distance distance metrics
#' @param method clustering method
#' @param cutoff number of clusters cutoff
#' @return cluster membership of the probes
ClusterMethylation <- function(methylation, distance='euclidean', method='ward.D2', cutoff=20) {
  cutoff = ifelse(cutoff < dim(methylation)[1], cutoff, dim(methylation)[1])
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
RegressMethylationOnRNAseq <- function(methylation, rnaseq, test.idx, method='glmnet', alpha=1) {
  stopifnot(length(methylation) == dim(rnaseq)[2])
  rnaseq.train <- t(data.matrix(rnaseq[, -test.idx]))
  rnaseq.test <- t(data.matrix(rnaseq[, test.idx]))
  methylation.train <- as.numeric(methylation[-test.idx])
  methylation.test <- as.numeric(methylation[, test.idx])
  if (method =='glmnet') {
      model.fit <- cv.glmnet(rnaseq.train, methylation.train, family='gaussian', 
                       alpha=alpha, type.measure='mse')
  } else if (method == 'vbsr') {
      model.fit <- vbsr(y=methylation.train, X=rnaseq.train, family='normal')
  } else {
      stop("method must be either glmnet or vbsr")
  }
  return(model.fit)
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

#' Extract coefficients that are significant from a vbsr object 
#' 
#' @param vbsr.fit a vbsr object
#' @return a data matrix of nonzero coefficient with row names as gene symbols
ExtractVBSRSignificantCoef <- function(vbsr.fit, m, coef.names) {
    p.val <- vbsr.fit$pval
    significant.idx <- qvalue(p.val, fdr.level=0.1)$significant
    coef.names <- coef.names[significant.idx]
    coef.significant <- vbsr.fit$beta[significant.idx]
    return(data.frame(gene=coef.names, coef=coef.significant))
}

#' Module network visualization of the methylation cluster, averge methylation profile
#' and predictive gene expression level
#' 
#' @param methylation a data matrix of methylation level (each probe is a column)
#' @param avemethyl a column vector of average methylation level across patient
#' @param rnaseq a data matrix of predictive gene RNAseq level (each probe is a column)
#' @param filename the path to which the image file (eps) will be saved
#' @return None (side effect: image saved at filename)
ModuleHeatmap <- function(methylation, avemethyl, rnaseq, filename=NULL, center=F, scale=F, coef=NULL) {
  rownames(methylation) <- NULL
  colnames(methylation) <- NULL
  rownames(rnaseq) <- NULL
  names(avemethyl) <- NULL
  ht1 <- Heatmap(methylation, col=inferno(10), width=unit(2.5, 'inches'), name='methylation', cluster_row=T, 
                 column_title='methylation')
  ht2 <- Heatmap(avemethyl, col=inferno(10), width=unit(0.1, 'inches'), name='average methylation', column_title="")
  if (center) {
    GEP <- scale(rnaseq, center=center, scale=scale)
    if(!is.null(coef)) GEP <- t(t(GEP) * coef)
  } else {
    GEP <- rnaseq 
  }
  ht3 <- Heatmap(GEP, col=inferno(10), width=unit(0.4, 'inches'),  name='GEP', column_title='GEP', cluster_columns=T)
  if(!is.null(filename)) {
    pdf(file=filename, width=6, height=10)
      draw(ht1 + ht2 + ht3)
    dev.off()
    print(paste('Image saved as', filename))
  } else {
    draw(ht1+ht2+ht3)
  }
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
                     cutoff=20, percent.test=0.3, alpha=1, seed=1000, center=F, scale=F, subset.method='glmnet', sample.code='01') {
  set.seed(seed)
  rnaseq.processed <- PreprocessRNAseq(rnaseq, sample.code=sample.code)
  rownames(rnaseq.processed$rnaseq) <- rnaseq.processed$g.geneid[,2]
  methylation.processed <- PreprocessMethylation(methylation, convert2M=convert2M,
                                                 subsetProbes=subsetProbes, sample.code=sample.code)
  processed <- matchColumns(rnaseq.processed$rnaseq, 
                            methylation.processed$methylation)
  methylation <- processed$methylation
  rnaseq <- processed$rnaseq
  methylprobes.geneid <- methylation.processed$m.geneid[,2]
  methylation <- AverageMethylation(methylation, methylprobes.geneid) # agregate the methylation level of probes of the same gene
  methylation <- data.matrix(methylation[,-1])
  rownames(methylation) <- methylation.processed$m.geneid[!duplicated(methylprobes.geneid), 1]
  
  test.idx <- sample(seq(dim(methylation)[2]), round(dim(methylation)[2]*percent.test))
  clustering <- ClusterMethylation(methylation, distance=distance, method=method, cutoff=cutoff)
  cluster <- clustering$cluster
  avemethyl <- AverageMethylation(methylation, cluster)
  for (i in seq(cutoff)) {
    if (sum(cluster==i) == 1) {
       print(paste('Only 1 probe ends up in cluster', i)) 
       print('Skipping to the next cluster')
       next
    }
    base.name <- paste0('cluster', i)
    print(paste('Running model for cluster', i))
    path.model <- paste0(datafolder, '/', base.name, '.RData')
    path.fig <- paste0(imagefolder, '/', base.name, '.eps')
    tryCatch({
    model <- RegressMethylationOnRNAseq(avemethyl[i, -1], rnaseq, test.idx, method=subset.method)
    if (subset.method=='glmnet') {
        coef <- ExtractNonzeroCoef(model)
    } else {
        coef <- ExtractVBSRSignificantCoef(model, m=dim(rnaseq)[1], coef.names=unlist(rnaseq.processed$g.geneid[,2]))
    }
    ordering <- order(abs(coef[-1,2]), decreasing=T)
    ordering.genes <- as.character(coef[-1,1][ordering])
    ModuleHeatmap(t(data.matrix(methylation[cluster==i, -test.idx])), 
                  t(data.matrix(avemethyl[i, -1][,-test.idx])), 
                  t(rnaseq[ordering.genes, -test.idx]), 
                  path.fig, center=center, scale=scale)
    save(model, coef, base.name, file=path.model)
    print(paste('Model saved at', path.model))
    }, error=function(e) {cat("ERROR :",conditionMessage(e), "\n")})
  }
  save(test.idx, clustering, avemethyl, 
       file=paste0(datafolder, '/', 'clusterInfo.RData'))
  print(paste('Cluster information saved at', path.model))
}
