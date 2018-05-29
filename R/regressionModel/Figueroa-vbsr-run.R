# File name: vbsr-run.R
# Author: Mamie Wang
# Date created: 04/23/2018
# Date modified: 04/23/2018

library(methods)
source('regressionModel.R')
source('../correlation_HM450/computeCorrelation.R')

args = commandArgs(trailingOnly=T)
if (length(args) < 3) {
  stop('At least five arguments must be supplied.', call.=F)
}

methylation = args[1]
rnaseq = args[2]
resultDir = args[3]
basename.dir <- basename(resultDir)

readTsv <- function(file.path) {
  data.raw <- read.csv(file.path, sep='\t', stringsAsFactor=T, header=F)
  data <- data.matrix(data.raw)[-1, -1]
  colnames(data) <- data.raw[1, -dim(data.raw)[2]]
  rownames(data) <- data.raw[-1, 1]
  return(data)
}

center <- T
scale <- T
dir.create(resultDir, showWarnings=F)
imagefolder <- paste0(basename.dir, '-plots/')
dir.create(file.path(resultDir, imagefolder), showWarnings=F)
imagefolder <- paste0(resultDir, imagefolder)
datafolder <- paste0(basename.dir, '-data/')
dir.create(file.path(resultDir, datafolder), showWarnings=F)
datafolder <- paste0(resultDir, datafolder)

print('Running regression model')
start.time <- proc.time()
RunModel2(methylation, RNAseq, imagefolder, datafolder)
print(proc.time() - start.time)
