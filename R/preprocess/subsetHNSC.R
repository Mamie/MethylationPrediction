# File name: vbsr-run.R
# Author: Mamie Wang
# Date created: 04/23/2018
# Date modified: 04/23/2018

library(methods)

args = commandArgs(trailingOnly=T)
if (length(args) < 2) {
  stop('At least two arguments must be supplied.', call.=F)
}

methylation = args[1]
subset = as.numeric(args[2])

load(methylation)

HM450.subset = data.HM450[1:subset,]
output = paste0(methylation, 'subset')
save(HM450.subset, file=output)
