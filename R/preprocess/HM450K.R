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
library(impute)
library(gdata)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop('At least two arguments must be supplied', call.=F)
}

in.file = args[1]
out.path = args[2]
subset.figueroa = args[3]

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
rm(columnnames, num.columns, redundant.idx)

print('The orginal data size is')
print(dim(data.HM450))
print('Filtering for methylation probes that are NA, chrX or ChrY')
data.HM450 <- data.HM450 %>%
  filter(!(Chromosome %in% c("X", "Y", NA)))
print('After removing chromosome X, Y, NA, the size is')
print(dim(data.HM450))

missingness <- apply(data.HM450[,-seq(4)], 1, function(x) mean(is.na(x)))
print(paste0('The number of probes with missingness >= 0.5 is ', sum(missingness >= 0.5)))
data.HM450 <- data.HM450[missingness < 0.5,]
print('After removing missing entries, the data size is')
print(dim(data.HM450))
data.HM450.imputed <- impute.knn(data.matrix(data.HM450[,-seq(4)]), k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
HM450.processed <- data.HM450.imputed$data
rownames(HM450.processed) <- c(data.HM450$Gene.Symbol)

print('converting to M value...')
beta2M <- function(beta) {
    log(beta/(1 - beta), base=2)
}
Beta2M <- Vectorize(beta2M)
HM450.processed <- apply(HM450.processed, 2, Beta2M)

print('Matching patients present in RNAseq...')
rnaseq <- read.csv('/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/LAML/processed/RNAseq-imputed.tsv', sep='\t', header=F, stringsAsFactor=F)
patients.rnaseq <- rnaseq[1,-dim(rnaseq)[2]]
patients.HM450 <- sapply(colnames(HM450.processed), function(x) toupper(substr(x, 1, 15)))
patients.matched <- which(patients.HM450 %in% colnames(rnaseq))
HM450.processed <- HM450.processed[, patients.matched]
print(dim(HM450.processed))

if(as.logical(subset.figueroa)) {
  n <- 16
  clusters <- list()
  for (i in seq(n)) {
    clusters[[i]] <- read.xls('/oak/stanford/groups/andrewg/users/szmamie/repos/MethylationPrediction/data/mmc3.xls', sheet=i, skip=1, header=T, stringsAsFactor=F)
  }
  cluster.probes <- list()
  for (i in seq(n)) {
    probe.names <- unique(clusters[[i]]$GENE.SYMBOL)
    cluster.probes[[i]] <- probe.names[probe.names!='']
  }
  
  findIdx <- function(probe.name) which(grepl(paste0(probe.name, '$|', probe.name, ';'), rownames(HM450.processed), fixed=F))
  FindIdx <- Vectorize(findIdx) 
  cluster.probes.all <- unique(unlist(cluster.probes))
  idx <- FindIdx(cluster.probes.all)
  probe.idx <- c(unlist(idx))
  print(paste0('The number of probes matching Figueroa paper are ', length(probe.idx)))
  HM450.processed <- HM450.processed[probe.idx,]
  print(dim(HM450.processed))
}

write.table(HM450.processed, file=out.path, sep='\t', quote=F, col.names=T, row.names=T)
