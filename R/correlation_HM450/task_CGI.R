# File name: task_CGI.R
# Author: Mamie Wang
# Date created: 03/01/2018
# Date modified: 03/01/2018
# Input: processed RNAseq and methylation data
# Output: correlation coefficients between matched pairs of RNAseq and methylation
# probes on the subset of promoter + CGI/shore region probes 
# (defined in notes/promoter-region-methylation-RNAseq-correlation.ipynb)

source('computeCorrelation.R')
load('../../data/processed/RNAseq.RData')
load('../../data/processed/HM450.RData')
CGI.probes <- read.csv('../../data/promoterRegionProbes/CGI.tsv', 
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(CGI.probes[,1])
correlation <- ComputeCorrelation(RNAseq, data.HM450, convert2M=T, subsetProbes=probe.id)
save(correlation, file='../../data/processed/CorrelationCGIM.RData')
correlation <- ComputeCorrelation(RNAseq, data.HM450, convert2M=F, subsetProbes=probe.id)
save(correlation, file='../../data/processed/CorrelationCGIB.RData')
