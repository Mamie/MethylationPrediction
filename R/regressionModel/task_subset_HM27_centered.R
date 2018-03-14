#!/usr/bin/env Rscript
# File name: task.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/13/2018

library(methods)
source('regressionModel.R')

load('../../data/processed/HM27.RData')
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
imagefolder <- '../../data/180313/HM27/moduleNetworkFigMcentersubset'
datafolder <- '../../data/180313/HM27/resultsMcentersubset'

RunModel(data.HM27, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale)

convert2M <- F
imagefolder <- '../../data/180313/HM27/moduleNetworkFigBcentersubset'
datafolder <- '../../data/180313/HM27/resultsBcentersubset'

RunModel(data.HM27, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale)


