library(methods)
source('../../notes/helperFunctions.R')

args = commandArgs(trailingOnly=T)
if (length(args) < 1) {
  stop('At least one arguments must be supplied.')
}

dataFolder = args[1]

WritePredictorInfo(dataFolder, 20)
