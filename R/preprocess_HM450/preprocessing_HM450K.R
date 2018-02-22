# File name: preprocessing_450K.R
# Author: Mamie Wang
# Date created: 02/22/2018
# Date last modified: 02/22/2018
# R version: 3.4.2
# Preprocess 450K methylation:
#   - remove NA, chrX or chrY

library(dplyr)

data.HM450 <- read.csv('../../data/methylationHM450/LAML.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt', 
                      sep='\t', header=T, stringsAsFactors=F)
num.columns <- dim(data.HM450)[2]
redundant.idx <- c(seq(7, num.columns, by=4), seq(8, num.columns, by=4), 
                   seq(9, num.columns, by=4))
data.HM450 <- data.HM450[-1,-redundant.idx]
data.HM450[,seq(5)] <- data.HM450[,c(1,3,4,5,2)]
columnnames <- colnames(data.HM450)
columnnames[seq(4)] <- c('ID', 'Gene.Symbol', 'Chromosome', 'Genomic.Coordinate')
colnames(data.HM450) <- columnnames
methProbes <- data.HM450[, seq(2)]
save(methProbes, file='../../data/processed/methylationProbes.RData')
rm(columnnames, num.columns, redundant.idx)
data.HM450 <- data.HM450 %>%
  filter(!(Chromosome %in% c("X", "Y", NA)))
missing.rm <- apply(data.HM450[seq(5, dim(data.HM450)[2])], 1, function(x) sum(is.na(x)) > 0)
data.HM450 <- data.HM450[!missing.rm,]
save(data.HM450, file='../../data/processed/HM450.RData')

