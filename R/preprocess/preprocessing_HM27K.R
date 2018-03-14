# File name: preprocessing_27K.R
# Author: Mamie Wang
# Date created: 03/07/2018
# Date last modified: 03/07/2018
# R version: 3.4.2
# Preprocess 27K methylation:
#   - remove NA, chrX or chrY

library(dplyr)

data.HM27 <- read.csv('../../data/methylationHM27/LAML.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt', 
                      sep='\t', header=T, stringsAsFactors=F)
num.columns <- dim(data.HM27)[2]
redundant.idx <- c(seq(7, num.columns, by=4), seq(8, num.columns, by=4), 
                   seq(9, num.columns, by=4))
data.HM27 <- data.HM27[-1,-redundant.idx]
data.HM27[,seq(5)] <- data.HM27[,c(1,3,4,5,2)]
columnnames <- colnames(data.HM27)
columnnames[seq(4)] <- c('ID', 'Gene.Symbol', 'Chromosome', 'Genomic.Coordinate')
colnames(data.HM27) <- columnnames
methProbes <- data.HM27[, seq(2)]
save(methProbes, file='../../data/processed/methylation27Probes.RData')
rm(columnnames, num.columns, redundant.idx)
data.HM27 <- data.HM27 %>%
  filter(!(Chromosome %in% c("X", "Y", NA)))
missing.rm <- apply(data.HM27[seq(5, dim(data.HM27)[2])], 1, function(x) sum(is.na(x)) > 0)
data.HM27 <- data.HM27[!missing.rm,]
save(data.HM27, file='../../data/processed/HM27.RData')

