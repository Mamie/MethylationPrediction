#!/usr/bin/env Rscript
# File name: task.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/13/2018

library(methods)
source('regressionModel.R')
source('../correlation_HM450/computeCorrelation.R')

load('../../data/processed/HM450.RData')
load('../../data/processed/RNAseq.RData')

CGI.probes <- read.csv('../../data/promoterRegionProbes/CGICorFdr10.tsv',
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(CGI.probes[,1])
subsetProbes <- probe.id

center <- T
scale <- T
convert2M <- T
percent.test <- 0.3
cutoff <- 20
alpha <- 1
distance <- 'euclidean'
method <- 'ward.D2'
imagefolder <- '../../data/180313/HM450/moduleNetworkFigMcenterScale'
datafolder <- '../../data/180313/HM450/resultsMcenterScale'

RunModel(data.HM450, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale)

convert2M <- F
imagefolder <- '../../data/180313/HM450/moduleNetworkFigBcenterScale'
datafolder <- '../../data/180313/HM450/resultsBcenterScale'

RunModel(data.HM450, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale)


