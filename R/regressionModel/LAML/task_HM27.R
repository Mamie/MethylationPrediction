#!/usr/bin/env Rscript
# File name: simpleTask.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/07/2018

library(methods)
source('regressionModel.R')
source('../correlation_HM450/computeCorrelation.R')

load('../../data/processed/HM27.RData')
load('../../data/processed/RNAseq.RData')
promoter.probes <- read.csv('../../data/promoterRegionProbes/TSS1500.tsv',
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(promoter.probes[,1])
subsetProbes <- probe.id

convert2M <- T
percent.test <- 0.3
cutoff <- 20
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

