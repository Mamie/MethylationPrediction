#!/usr/bin/env Rscript
# File name: task.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/07/2018

library(methods)
source('regressionModel.R')

load('../../data/processed/HM450.RData')
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
imagefolder <- '../../data/180307_M/HM450/moduleNetworkFigM'
datafolder <- '../../data/180307_M/HM450/resultsM'

RunModel(data.HM450, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha)

convert2M <- F
imagefolder <- '../../data/180307_M/HM450/moduleNetworkFigB'
datafolder <- '../../data/180307_M/HM450/resultsB'

RunModel(data.HM450, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha)


