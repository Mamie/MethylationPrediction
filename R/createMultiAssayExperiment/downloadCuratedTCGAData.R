# file name: downloadCuratedTCGAData.R
# Author: Mamie Wang
# Date created: 03/13/2018
# Date last modified: 03/13/2018
# Download curated TCGA HNSC data and create a MultiAssayExperiment Object from
# Bioconductor package curatedTCGAData

library(curatedTCGAData)
curatedTCGAData(diseaseCode = c("HNSC"), 
		assays = c("Methylation", "RNASeq2GeneNorm"),
                dry.run=F)
curatedTCGAData(diseaseCode = c("SKCM"), 
                assays = c("Methylation", "RNASeq2GeneNorm"),
                dry.run=F)
