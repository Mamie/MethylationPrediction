# File name: task_Mvalue.R
# Author: Mamie Wang
# Date created: 02/28/2018
# Date modified: 02/28/2018
# Input: processed RNAseq and methylation data
# Output: correlation coefficients between matched pairs of RNAseq and methylation probes

source('computeCorrelation.R')
load('../../data/processed/RNAseq.RData')
load('../../data/processed/HM450.RData')

correlation <- ComputeCorrelation(RNAseq, data.HM450, T)
save(correlation, file='../../data/processed/CorrelationMvalue.RData')
