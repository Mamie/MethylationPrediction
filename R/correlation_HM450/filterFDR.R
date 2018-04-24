library(qvalue)
source('computeCorrelation.R')

args = commandArgs(trailingOnly=T)

if (length(args) < 4) {
   stop('At least four arguments must be supplied.')
}

methylation = args[1]
rnaseq = args[2]
outpath = args[3]
fdr = as.numeric(args[4])

print('loading methylation data...')
load(methylation)
print('loading rnaseq data...')
load(rnaseq)

CGI.probes <- read.csv("/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/promoterRegionProbes/CGI.tsv",
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(CGI.probes[,1])

print('computing correlation between methylation and gene expression')
correlation <- ComputeCorrelation(RNAseq, data.HM450, convert2M=T, subsetProbes=probe.id)

negcor <- correlation$correlation < 0
p.value <- correlation$p.value
q.obj <- qvalue(p = p.value, fdr.level=fdr)
probes.significant <- matrix(rownames(correlation)[q.obj$significant & negcor], ncol=1)
colnames(probes.significant) <- c('probe')
print(paste0('writing the results to ', outpath))
write.table(probes.significant, file=outpath,
   sep='\t', quote=F, row.names=F)
print(paste0('The number of probes that are in the CGI region and are negatively correlated with GEP with fdr < 0.1 is ', dim(probes.significant)[1]))
