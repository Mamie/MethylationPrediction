# Author: Mamie Wang
# Date created: 02/05/2018
# Date last modified: 02/06/2018
# R version: 3.4.2
# Run lasso model for each methylation cluster using mRNA-seq genes as predictors
#   - group methylation beta-value within each cluster
#   - filter mRNA-seq for probes that corresponds to genes with keyword (methyl,
#   transcription)

library(glmnet)
library(dplyr)

load('../data/processed/clusterTCGA.RData')
load('../data/processed/MethylationFiltered.RData')
load('../data/processed/filteredGEP.RData')
load('../data/processed/RNAseq.RData')
load('../data/processed/transcriptionRelatedGenes.RData')

# Group beta values within clusters
print('Average beta values within clusters...')
probeCluster <- genes %>%
  dplyr::select(probe, cluster) %>%
  dplyr::distinct()
M.grouped <- data.frame(M.filtered) %>%
  dplyr::mutate(probe = as.character(rownames(M.filtered))) %>%
  dplyr::left_join(probeCluster) %>%
  dplyr::group_by(cluster) %>%
  dplyr::select(-probe) %>%
  dplyr::summarise_all(funs(mean)) %>%
  na.omit()
save(M.grouped, file='../data/processed/averageMethylationProfile.RData')

i = 1
methClusterNums <- unique(genes$cluster)
for (i in seq(methClusterNums)) {
  M.cluster <- data.frame(M.filtered) %>%
    dplyr::mutate(probe = as.character(rownames(M.filtered))) %>%
    dplyr::left_join(probeCluster) %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(cluster==i)
  save(M.cluster, file=paste0('../data/processed/Mcluster/Mcluster', i, '.RData'))
}

# Filtering for RNAseq genes with 'methyl' or 'transcription' keyword
print('Obtaining RNAseq genes with specific keywords')
filtered.genes <- as.character(terms.filtered$name)
RNAseq.genes <- sapply(as.character(G.filtered[,1]), function(x) strsplit(x, '|', fixed=T)[[1]][1])
G <- G.filtered[RNAseq.genes %in% filtered.genes,] # 1629  166
row.names(G) <- RNAseq.genes[RNAseq.genes %in% filtered.genes]
G <- G[,-1]
dim(G) # 1629  166
genes.GEP <- RNAseq.genes[RNAseq.genes %in% filtered.genes]


# Run lasso for each cluster
print('Running lasso for each cluster...')
for (idx in seq(methClusterNums)) {
  cluster <- genes %>%
    filter(cluster==idx)
  cluster.probes <- cluster %>%
    dplyr::select(probe) %>% dplyr::distinct() %>% unlist()
  excludeGenes <- cluster %>%
    dplyr::select(Name) %>% unlist()
  print(paste0('Found overlapping genes between mRNA-seq and methylation array ',
               length(excludeGenes), ' genes'))
  
  G.cluster <- G[!(genes.GEP %in% excludeGenes),]
  #G.cluster <- G[genes.GEP,]
  print('Dimension of the cluster data:')
  print(dim(G.cluster))
  
  # 70 % / 30 % training / test split
  set.seed(100)
  patient.test <- sample(seq(dim(G.cluster)[2]), round((dim(G.cluster)[2] - 1)*0.3))
  G.clustertrain <- t(data.matrix(G.cluster[, -patient.test]))
  G.clustertest <- t(data.matrix(G.cluster[,patient.test]))
  M.cluster <- M.grouped[idx,][,-1]
  M.clustertrain <- as.numeric(M.cluster[, -patient.test])
  M.clustertest <- as.numeric(M.cluster[, patient.test])
  
  # Lasso model
  for (i in 0:10) {
    assign(paste('fit.cluster', i, sep=''), cv.glmnet(G.clustertrain, M.clustertrain, 
                                                      family='gaussian', alpha=i/10, 
                                                      type.measure='mse'))
  }
  
  # nonzero-coefficients
  coef.lasso <- data.matrix(predict(fit.cluster10, type='coef', s='lambda.min'))
  sum(coef.lasso != 0) # 41 
  nonzero.coefraw <- rownames(coef.lasso)[coef.lasso != 0][-1]
  nonzero.coef <- sapply(nonzero.coefraw, function(x) strsplit(x, '|', fixed=T)[[1]][1])
  names(nonzero.coef) <- NULL
  
  # coefficients
  nonzero.coefvalues <- coef.lasso[coef.lasso != 0][-1]
  names(nonzero.coefvalues) <- nonzero.coef
  nonzero.coefvalues <- sort(nonzero.coefvalues)
  nonzero.coefannot <- data.frame(name=names(nonzero.coefvalues), value=nonzero.coefvalues)  %>%
    left_join(terms.filtered)
  save(nonzero.coefannot, file=paste0('../data/processed/coefs/coef', idx, '.RData'))
  G.clusternonzero <- G.cluster[rownames(G.cluster) %in% nonzero.coefraw, ]
  save(G.clusternonzero, file=paste0('../data/processed/GEP/GEP', idx, '.RData'))
}
