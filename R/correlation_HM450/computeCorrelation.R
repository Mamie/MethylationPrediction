#!/usr/bin/env Rscript
# File name: computeCorrelation.R
# Author: Mamie Wang
# Date created: 02/27/2018
# Date modified: 02/27/2018
# Input: 
#  - a hash map from methylation probe id to gene probe id
#  - methylation data matrix with row names as methylation probe id
#  - RNAseq data matrix with row names as RNAseq probe id and same number of columns
#  - of the methylation data matrix
# Output: a vector of correlation coefficients for each pair of methylation and RNAseq

MethylationGEPCorrelation <- function(hash.map, methylation, rnaseq) {
  stopifnot(dim(methylation)[2] == dim(rnaseq)[2])
  methylation.id <- row.names(methylation)
  rnaseq.id <- row.names(rnaseq)
  rnaseq.idmatched <- unlist(sapply(methylation.id, function(x) hash.map[[x]]))
  rnaseq.index <- unlist(sapply(rnaseq.idmatched, function(x) { 
    match <- rnaseq.id == x
    if (sum(match)== 0) return(NA)
    return(which(rnaseq.id==x))
    }))
  pairs.index <- cbind(seq(length(methylation.id)), rnaseq.index)[which(!is.na(rnaseq.index)),]
  
  unlist(apply(pairs.index, 1, function(x) cor(methylation[x[1],], rnaseq[x[2],])))
}
