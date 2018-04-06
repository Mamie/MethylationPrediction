#!/usr/bin/env Rscript
# File name: task.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date modified: 03/14/2018

library(methods)
source('regressionModel.R')

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop('At least two arguments must be supplied', call.=F)
}
data.dir = args[1]

load(paste0(data.dir, '/HM450.RData'))
load(paste0(data.dir, '/RNAseq.RData'))

CGI.probes <- read.csv(paste0(data.dir, '/CGICorFdr10.tsv'),
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
imagefolder <- paste0(data.dir, '/moduleNetworkFigMcenterScale')
datafolder <- paste0(data.dir, '/resultsMcenterScale')

RunModel(data.HM450, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale)

convert2M <- F
imagefolder <- paste0(data.dir, '/moduleNetworkFigBcenterScale')
datafolder <- paste0(data.dir, '/resultsBcenterScale')

RunModel(data.HM450, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale)


