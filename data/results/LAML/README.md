# Methylation clustering and GSEA analysis results

Beta values of the HM27 probes with the same gene name are averaged and hierarchical clustering 
were performed on the resulting data matrix. 20 clusters were selected as cutoff criterion. Probe ID, corresponding
HUGO gene symbol and description for each of the cluster is located in `vbsr-M-plots` named as `cluster$i.csv`. 

The probes in each cluster were investigated by computing the overlaps with other gene sets in MSigDB (C2, C3), 
the results were located in `vbsr-M-plots` named as `cluster$i.xls`.
(Cluster 2, 6, 10, 12, 18, 20 has no overlaps with C2 and C3 gene sets.)

# VBSR subset selection results

VBSR were used to select subset of predictors. Gene symbol, beta coefficients and HUGO gene symbol
description of the selected predictors are located `vbsr-M-data` named as `cluster$iPredictor.csv`

The module network plot of the methylation cluster and predictive genes can be found in `vbsr-M-plots/vbsrMHM27.pdf`
