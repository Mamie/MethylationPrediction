# File name: 020618_correlation.R
# Author: Mamie Wang
# Date created: 02/06/2018
# Date last modified: 02/06/2018
# R version: 3.4.2
# 1. Correlation between mRNA-seq level and methylation beta value
# 2. Compare clustering between consensus clustering and hclust k = 20

library(dplyr)
library(tidyr)
library(fpc)
library(ConsensusClusterPlus)

# Correlation between mRNA-seq and methylation level
load("../data/processed/Beta190.RData")
load("../data/processed/RNAseq.RData")
load("../data/processed/methylationProbes.RData")

patients.intersect <- intersect(colnames(RNAseq), colnames(Beta190))
Beta190 <- Beta190[, patients.intersect]
RNAseq <- RNAseq[, patients.intersect]

geneSymbols <- strsplit(methProbes$Gene.Symbol, ";")
max_len <- max(sapply(geneSymbols, length))
for(i in 1:max_len){
  methProbes[, paste0("Gene.Symbol", i)] <- sapply(geneSymbols, "[", i)
}

# get methylation probes - gene correspondence
genes <- methProbes %>%
  dplyr::select(-Gene.Symbol) %>%
  tidyr::gather(Gene.Symbol, Name, Gene.Symbol1:Gene.Symbol22) %>%
  dplyr::filter(!is.na(Name)) 

# get mRNA-seq genes correspondence
ExtractGeneSymbol <- function(x) strsplit(x, '|', fixed=T)[[1]][1]
genes.RNAseq <- data.frame(rownum=seq(dim(RNAseq)[1]), 
                           geneid=sapply(RNAseq$gene, 
                                         ExtractGeneSymbol)) %>%
  mutate_if(is.factor, as.character)
  
# join the two table and filter for pairs of methylation and gene probes
meth.gene.Pair <- genes %>%
  left_join(genes.RNAseq, by=c('Name'='geneid')) %>%
  filter(!is.na(rownum)) %>%
  select(ID, rownum, Name)

computeCorrelation <- function(x) {
  cor(unlist(Beta190[x[1], -1]), unlist(RNAseq[as.numeric(x[2]), -1]))
}
correlation <- apply(meth.gene.Pair, 1, computeCorrelation)
save(correlation, file='../data/processed/correlation.RData')
hist(correlation)

# Compare clustering between hclust k = 20 and consensus clustering

load('../data/processed/MethylationFiltered.RData')
distances <- dist(M.filtered) # euclidean distance
hclust.M <- hclust(distances, method='ward.D2')
fit.hclust <- cutree(hclust.M, 20)

M.filtered <- sweep(M.filtered, 1, apply(M.filtered, 1, median, na.rm=T))
pItem <- 0.80
pFeature <- 0.80
maxK <- 20
reps <- 50
clustAlg <- 'hc'
distance <- 'pearson'
title <- '../data/figures/Methylation'
seed <- 10000
results <- ConsensusClusterPlus(d=t(M.filtered), maxK=maxK, reps=reps, pItem=pItem, 
                               pFeature=pFeature, title=title, clusterAlg=clustAlg,
                               distance=distance, seed=seed, plot='png') #'pngBMP' for large datasets
fit.consensusclust <- results[[20]]$consensusClass
save(distances, hclust.M, fit.consensusclust, file='../data/processed/cluster.RData')
cluster.stats(distances, fit.hclust, fit.consensusclust)
