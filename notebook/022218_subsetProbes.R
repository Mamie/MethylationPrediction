# filename: 022218_subsetProbes.R
# author: Mamie Wang
# date created: Feb 22 2018
# date modified: Feb 22 2018
# This script subsets for the promoter region methylation probes
# And compute the correlation with corresponding genes
library(dplyr)

annotations <- read.csv('../data/HumanMethylation450_15017482_v1-2.csv', 
                        sep=',', stringsAsFactors=F, skip=7, header=T)
row.names(annotations) <- annotations$IlmnID
promoter.region <- annotations$Phantom != ''
probes.promoter <- rownames(annotations)[promoter.region]

load('../data/processed/methylationMvalue.RData')
probes.promoter <- intersect(probes.promoter, rownames(methylation.M)) # 178

methylation.promoters <- methylation.M[rownames(methylation.M) %in% probes.promoter,]

load("../data/processed/RNAseq.RData")
load("../data/processed/methylationProbes.RData")
patients.intersect <- intersect(colnames(RNAseq), colnames(methylation.promoters))
methylation.promoters <- methylation.promoters[, patients.intersect]
RNAseq.genes <- RNAseq$gene
RNAseq <- RNAseq[, patients.intersect]

methProbes <- methProbes[methProbes$ID %in% probes.promoter,] # 2534
geneSymbols <- strsplit(methProbes$Gene.Symbol, ";")
max_len <- max(sapply(geneSymbols, length))
for(i in 1:max_len){
  methProbes[, paste0("Gene.Symbol", i)] <- sapply(geneSymbols, "[", i)
}

# get methylation probes - gene correspondence
genes <- methProbes %>%
  dplyr::select(-Gene.Symbol) %>%
  tidyr::gather(Gene.Symbol, Name, Gene.Symbol1:Gene.Symbol20) %>%
  dplyr::filter(!is.na(Name)) 

# get mRNA-seq genes correspondence
ExtractGeneSymbol <- function(x) strsplit(x, '|', fixed=T)[[1]][1]
genes.RNAseq <- data.frame(rownum=seq(dim(RNAseq)[1]), 
                           geneid=sapply(RNAseq.genes, 
                                         ExtractGeneSymbol)) %>%
  mutate_if(is.factor, as.character)

# join the two table and filter for pairs of methylation and gene probes
meth.gene.Pair <- genes %>%
  left_join(genes.RNAseq, by=c('Name'='geneid')) %>%
  filter(!is.na(rownum)) %>%
  select(ID, rownum, Name)

computeCorrelation <- function(x) {
  cor(unlist(methylation.promoters[x[1],]), unlist(RNAseq[as.numeric(x[2]),]))
}
correlation <- apply(meth.gene.Pair, 1, computeCorrelation)
save(correlation, file='../data/processed/correlation.RData')
hist(correlation, main='promoter region correlation')
