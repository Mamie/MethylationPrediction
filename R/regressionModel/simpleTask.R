#!/usr/bin/env Rscript
# File name: simpleTask.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/07/2018

library(methods)
source('regressionModel.R')

load('../../data/processed/subsetHM450.RData')
load('../../data/processed/subsetRNAseq.RData')

rnaseq <- data.rnaseq
methylation <- data.HM450
convert2M <- T
subsetProbes <- NULL
percent.test <- 0.3
cutoff <- 2
alpha <- 1
distance <- 'euclidean'
method <- 'ward.D2'
imagefolder <- '../../data/180307_M/moduleNetworkFig'
datafolder <- '../../data/180307_M/results'

RunModel(methylation, rnaseq, imagefolder, datafolder, convert2M=convert2M, 
         subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha)