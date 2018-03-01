# File name: task_promoter.R
# Author: Mamie Wang
# Date created: 03/01/2018
# Date modified: 03/01/2018
# Input: processed RNAseq and methylation data
# Output: correlation coefficients between matched pairs of RNAseq and methylation
# probes on the subset of promoter region probes 
# (defined in notes/promoter-region-methylation-RNAseq-correlation.ipynb)

source('computeCorrelation.R')
load('../../data/processed/RNAseq.RData')
load('../../data/processed/HM450.RData')
promoter.probes <- read.csv('../../data/promoterRegionProbes/TSS1500.tsv', 
                            sep='\t', header=T, stringsAsFactors=F)
probe.id <- as.character(promoter.probes[,1])
correlation <- ComputeCorrelation(RNAseq, data.HM450, T, probe.id)
save(correlation, file='../../data/processed/CorrelationPromoterM.RData')
correlation <- ComputeCorrelation(RNAseq, data.HM450, F, probe.id)
save(correlation, file='../../data/processed/CorrelationPromoterB.RData')


