library(methods)
source('../../notes/helperFunctions.R')

args = commandArgs(trailingOnly=T)
if (length(args) < 2) {
  stop('At least two arguments must be supplied.')
}

dataFolder = args[1]
methylationProbes = args[2]

print(paste0('Loading ', methylationProbes))
load(methylationProbes)
clusterInfo <- paste0(dataFolder, '/', 'clusterInfo.RData')

print(paste0('Loading ', clusterInfo, '...'))
load(clusterInfo)

dir.create(file.path(dataFolder, 'methylationClusters'), showWarnings=F)
WriteClusterInfor(clustering$cluster, methProbes, file.path(dataFolder, 'methylationClusters'))
