# Author: Mamie Wang
# Date created: 02/05/2018
# Date last modified: 02/05/2018
# R version: 3.4.2
# Module network visualization of the clusters

library(ComplexHeatmap)
library(viridis)

ht_global_opt(heatmap_column_names_gp = gpar(fontsize = 6), 
              heatmap_legend_title_gp= gpar(fontsize = 8),
              heatmap_legend_labels_gp=gpar(fontsize = 6),
              heatmap_row_names_gp = gpar(fontsize = 6))


# average methylation level for each patient within each cluster
load("../data/processed/averageMethylationProfile.RData")

methClusterNums = 20
i = 1
for (i in seq(methClusterNums)) {
  print(paste('Plotting for cluster', i))
  # load Lasso model coefficients for cluster
  load(paste0("../data/processed/coefs/coef", i, ".RData"))
  
  # load GEP level for the Lasso model nonzero coefficients for cluster
  load(paste0("../data/processed/GEP/GEP", i, ".RData"))
  
  # load methylation level for cluster
  load(paste0("../data/processed/Mcluster/Mcluster", i, ".RData"))
  
  M <- t(data.matrix(M.cluster[, seq(dim(M.cluster)[2]-2)]))
  rownames(M) <- NULL
  colnames(M) <- NULL
  ht1 <- Heatmap(M, col=inferno(256), name='methylation', cluster_row=T, column_title='methylation')
  
  ha1 <- t(data.matrix(M.grouped[i, -1]))
  rownames(ha1) <- NULL
  colnames(ha1) <- NULL
  ht2 <- Heatmap(ha1, col=inferno(256), width=unit(0.5, 'cm'), name='Average methylation')
  
  GEP <- scale(t(G.clusternonzero[as.character(nonzero.coefannot[,1]), ]),
               center=TRUE, scale=FALSE)
  rownames(GEP) <- NULL
  colnames(GEP) <- NULL
  ht3 <- Heatmap(GEP, col=inferno(256), name='GEP', column_title='GEP', cluster_columns=F)
  
  print(paste0('../data/figures/moduleVisualization/cluster', i, '.png'))
  png(filename=paste0('../data/figures/moduleVisualization/cluster', i, '.png'))
  draw(ht1 + ht2 + ht3)
  dev.off()
}
