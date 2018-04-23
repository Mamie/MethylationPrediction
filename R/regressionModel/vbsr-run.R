# File name: vbsr-run.R
# Author: Mamie Wang
# Date created: 04/23/2018
# Date modified: 04/23/2018

library(methods)
source('regressionModel.R')
source('../correlation_HM450/computeCorrelation.R')

args = commandArgs(trailingOnly=T)
if (length(args) < 4) {
  stop('At least four arguments must be supplied.', call.=F)
}

cancer = args[1]
methylation = args[2]
rnaseq = args[3]
convert2M = as.numeric(args[4])

load(methylation)
load(rnaseq)

if(exists("data.HM27")) {
   methylation = data.HM27
   hm = 'hm27'
   print('Found data.HM27')
   rm(data.HM27)
} else if(exists("data.HM450")) {
   methylation = data.HM450
   hm = 'hm450'
   print('Found data.HM450')
   rm(data.HM450)
} else {
   methylation = HM450.subset
   hm = 'hm450'
#   stop("Methylation data not found")
}

CGI.probes <- read.csv('../../data/promoterRegionProbes/CGI.tsv',
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(CGI.probes[,1])
subsetProbes <- probe.id

center <- T
scale <- T
convert2M <- as.logical(convert2M)
percent.test <- 0.3
cutoff <- 20
alpha <- 1
subset.method <- 'vbsr'
distance <- 'euclidean'
method <- 'ward.D2'
dir.create(cancer, showWarnings=F)
imagefolder <- paste0('vbsr-', hm, ifelse(convert2M, '-M', '-B'), '-plots/')
dir.create(file.path(cancer, imagefolder), showWarnings=F)
imagefolder <- paste0(cancer, imagefolder)
datafolder <- paste0('vbsr-', hm, ifelse(convert2M, '-M', '-B'), '-data/')
dir.create(file.path(cancer, datafolder), showWarnings=F)
datafolder <- paste0(cancer, datafolder)

print('Running regression model')
start.time <- proc.time()
RunModel(methylation, RNAseq, imagefolder, datafolder, convert2M=convert2M, subsetProbes=subsetProbes, distance=distance, method=method, 
         cutoff=cutoff, percent.test=percent.test, alpha=alpha, center=center, scale=scale, subset.method=subset.method)
print(proc.time() - start.time)
