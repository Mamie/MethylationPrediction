#!/usr/bin/env Rscript
# File name: simpleTask.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/07/2018

library(methods)
source('regressionModel.R')

load('../../data/processed/HM27.RData')
load('../../data/processed/RNAseq.RData')

convert2M <- T
subsetProbes <- NULL
percent.test <- 0.3
cutoff <- 2
alpha <- 1
distance <- 'euclidean'
method <- 'ward.D2'
imagefolder <- '../../data/180307_M/HM27/moduleNetworkFigM'
datafolder <- '../../data/180307_M/HM27/resultsM'

RunModel(data.HM27, RNAseq, imagefolder, datafolder, convert2M=convert2M, 
         subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha)

convert2M <- F
imagefolder <- '../../data/180307_M/HM27/moduleNetworkFigB'
datafolder <- '../../data/180307_M/HM27/resultsB'
RunModel(data.HM27, RNAseq, imagefolder, datafolder, convert2M=convert2M, 
         subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha)

