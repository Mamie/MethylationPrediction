# File name: preprocessing_450K.R
# Author: Mamie Wang
# Date created: 02/22/2018
# Date last modified: 03/14/2018
# R version: 3.4.2
# Input: the HM450 data file, and output directory for preprocessed data
# Preprocess 450K methylation:
#   - remove NA, chrX or chrY

library(dplyr)
library(methods)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop('At least two arguments must be supplied', call.=F)
}

in.file = args[1]
out.dir = args[2]

print('Reading in the data set...')
data.HM450 <- read.csv(in.file, sep='\t', header=T, stringsAsFactors=F)
print('Preprocessing to remove redudant columns...')
num.columns <- dim(data.HM450)[2]
redundant.idx <- c(seq(7, num.columns, by=4), seq(8, num.columns, by=4), 
                   seq(9, num.columns, by=4))
data.HM450 <- data.HM450[-1,-redundant.idx]
data.HM450[,seq(5)] <- data.HM450[,c(1,3,4,5,2)]
columnnames <- colnames(data.HM450)
columnnames[seq(4)] <- c('ID', 'Gene.Symbol', 'Chromosome', 'Genomic.Coordinate')
colnames(data.HM450) <- columnnames
methProbes <- data.HM450[, seq(2)]
save(methProbes, file=paste0(out.dir, 'methylationProbes.RData'))
rm(columnnames, num.columns, redundant.idx)

print('The orginal data size is')
print(dim(data.HM450))
print('Filtering for methylation probes that are NA, chrX or ChrY')
data.HM450 <- data.HM450 %>%
  filter(!(Chromosome %in% c("X", "Y", NA)))
print('After removing chromosome X, Y, NA, the size is')
print(dim(data.HM450))
missing.rm <- apply(data.HM450[seq(5, dim(data.HM450)[2])], 1, function(x) sum(is.na(x)) > 0)
data.HM450 <- data.HM450[!missing.rm,]
print('After removing missing entries, the data size is')
print(dim(data.HM450))
save(data.HM450, file=paste0(out.dir, 'HM450.RData'))
