# File name: task_FDR10.R
# Author: Mamie Wang
# Date created: 03/13/2018
# Date modified: 03/13/2018
# Input: processed RNAseq and methylation data
# Output: methylation probes with correlation coefficients between matched pairs of RNAseq and methylation
# probes on the subset of promoter + CGI/shore region probes less than FDR 0.1

library(qvalue)
source('computeCorrelation.R')

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop('At least one arguments must be supplied', call.=F)
}
data.dir = args[1]

load(paste0(data.dir, '/RNAseq.RData'))
load(paste0(data.dir, '/HM450.RData'))

q <- 0.1

CGI.probes <- read.csv('../../data/promoterRegionProbes/CGI.tsv', 
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(CGI.probes[,1])
correlation <- ComputeCorrelation(RNAseq, data.HM450, convert2M=T, subsetProbes=probe.id)
p.value <- correlation$p.value
q.obj <- qvalue(p = p.value, fdr.level=q)
probes.significant <- matrix(rownames(correlation)[q.obj$significant], ncol=1)
colnames(probes.significant) <- c('probe')
write.table(probes.significant, file=paste0(data.dir, '/CGICorFdr10.tsv'), sep='\t', quote=F, row.names=F)
print(dim(probes.significant)[1])

