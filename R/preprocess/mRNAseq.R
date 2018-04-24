# file name: mRNAseq.R
# Author: Mamie Wang
# Date created: 03/14/2018
# Date last modified: 04/23/2018
# Input: mRNAseq data file, output directory for processed data

library(methods)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop('At least two arguments must be supplied')
}
in.file = args[1]
out.dir = args[2]

RNAseq <- read.csv(in.file, sep='\t', header=T, stringsAsFactors=F)
print('Original data size is')
print(dim(RNAseq))
print('Preprocessing mRNA-seq data...')
RNAseq <- na.omit(RNAseq)
print('After removing probes with NA, the data size is')
print(dim(RNAseq))
FormatPatientID <- function(x) tolower(x)
colnames(RNAseq) <- sapply(colnames(RNAseq), FormatPatientID)
save(RNAseq, file=out.dir)
