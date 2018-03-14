#!/usr/bin/env Rscript
# File name: computeCorrelation.R
# Author: Mamie Wang
# Date created: 02/27/2018
# Date modified: 03/6/2018

library(hashmap)

#' Construct a vector of correlation coefficients for each pair of methylation and RNAseq.
#' 
#' @param hash.map a hash map from methylation probe id to gene probe id
#' @param methylation methylation data matrix with row names as methylation probe id
#' @param rnaseq RNAseq data matrix with row names as RNAseq probe id and 
#' same number of columns of the methylation data matrix
#' @param m.geneid A data frame with first column methylation probe id and second column
#' corresponding gene symbol
#' @return a data frame with two columns:
#'  - correlation coefficients all possible matching pair of methylation and RNAseq
#'  - the gene of the pair

MethylationGEPCorrelation <- function(hash.map, methylation, rnaseq, m.geneid) {
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
  colnames(pairs.index) <- NULL
  correlation <- apply(pairs.index, 1, function(x) {
			test <- cor.test(as.numeric(methylation[x[1],]), 
                                         as.numeric(rnaseq[x[2],]))
			return(c(test$statistic, test$p.value))
			})
  m2g.hashmap <- hashmap(m.geneid[,1], m.geneid[,2])                                                          
  gene <- unlist(sapply(colnames(correlation), function(x) m2g.hashmap[[x]]))
  return(data.frame(gene=gene, correlation=correlation[1,], p.value=correlation[2,]))
}

#' Construct a hash map from methylation probe id to RNAseq probe id
#' 
#' @param m.geneid a data matrix with first column methylation probe id and 
#' second column corresponding gene
#' @param g.geneid a data matrix with first column RNAseq probe id and 
#' second column corresponding gene (one-to-one correspondence)
#' @return a hash map with key methylation probe id and values RNAseq probe id
ConstructM2GHashMap <- function(m.geneid, g.geneid) {
  g.hashmap <- hashmap(g.geneid[,2], g.geneid[,1])
  hashmap(m.geneid[,1], sapply(m.geneid[,2], function(x) g.hashmap[[x]]))
}

#' Preprocess the RNAseq data for input of MethylationGEPCorrelation
#' 
#' @param rnaseq RNAseq data frame with first column gene id
#' @return a list with row number-gene id data matrix, and data matrix with only RNAseq level
PreprocessRNAseq <- function(rnaseq) {
  geneid <- sapply(rnaseq$gene, function(x) strsplit(x, '|', fixed=T)[[1]][1])
  notna <- geneid != '?'
  notduplicated <- !duplicated(geneid)
  filtered <- notna & notduplicated
  rnaseq.filtered <- rnaseq[filtered, -1]
  rownames(rnaseq.filtered) <- seq(sum(filtered))
  g.geneid <- cbind(as.character(seq(sum(filtered))), geneid[filtered])
  rownames(g.geneid) <- NULL
  return(list(g.geneid=g.geneid, rnaseq=rnaseq.filtered))
}

#' Preprocess the methylation data for input of MethylationGEPCorrelation
#' 
#' @param methylation Methylation data frame (first column probe id and second 
#' @param convert2M whether converting to m value
#' column gene id) 
#' @param subsetProbes a character vector of the subset of the probe ids
#' @return a list with probe id - gene id data matrix and data matrix with only
#' methylation level
PreprocessMethylation <- function(methylation, convert2M=F, subsetProbes=NULL) {
  probeid <- methylation$ID
  geneid <- methylation$Gene.Symbol
  if(!is.null(subsetProbes)) {
    insubset <- probeid %in% subsetProbes
  } else {
    insubset <- rep(T, length(probeid))
  }
  notna <- !is.na(geneid)
  filtered <- notna & insubset
  methylation.filtered <- data.matrix(methylation[filtered, -seq(4)])
  if (convert2M) {
    methylation.filtered <- Beta2M(methylation.filtered, 10)
  }
  colnames(methylation.filtered) <- sapply(colnames(methylation.filtered),
                                           function(x) tolower(substr(x, 1, 12)))
  rownames(methylation.filtered) <- probeid[filtered]
  m.geneid <- cbind(probeid[filtered], geneid[filtered])
  rownames(m.geneid) <- NULL
  return(list(m.geneid=m.geneid, methylation=methylation.filtered))
}

#' Match the columns between RNAseq and Methylation array
#' 
#' @param rnaseq RNAseq data matrix
#' @param methylation Methylation data matrix
#' @return a list containing RNAseq and methylation with matching columns
matchColumns <- function(rnaseq, methylation) {
  intersection <- intersect(colnames(rnaseq), colnames(methylation))
  intersection.rnaseq <- sapply(intersection, function(x) which(colnames(rnaseq) %in% x))
  intersection.methylation <- sapply(intersection, function(x) which(colnames(methylation) %in% x))
  list(rnaseq=rnaseq[,intersection.rnaseq], 
       methylation=methylation[,intersection.methylation])
}

#' Main function to run the correlation coefficient function
#'
#' @param rnaseq RNAseq data frame
#' @param methylation methylation data frame
#' @param convet2M whether converting to M value for methylation
#' @param subsetProbes a subset of probe id
#' @return a vector of correlation coefficient between matching pairs of rnaseq 
#' and methylation probes
ComputeCorrelation <- function(rnaseq,  methylation, convert2M=F, subsetProbes=NULL) {
  rnaseq.processed <- PreprocessRNAseq(rnaseq)
  methylation.processed <- PreprocessMethylation(methylation, convert2M=convert2M, subsetProbes=subsetProbes)
  hash.map <- ConstructM2GHashMap(methylation.processed$m.geneid, 
                      rnaseq.processed$g.geneid)
  processed <- matchColumns(rnaseq.processed$rnaseq, 
                            methylation.processed$methylation)
  correlation <- MethylationGEPCorrelation(hash.map, processed$methylation, 
                                           processed$rnaseq, 
                                           methylation.processed$m.geneid)
  return(correlation)
}

#' Convert to the methylation beta value to M value 
#' 
#' @param methylation data matrix of methylation beta value
#' @param threshold threshold (postive number) at which to truncate M values 
#' @return corresponding data matrix of M value
Beta2M <- function(methylation, threshold) {
  M <- apply(methylation, 2, function(beta) log(beta/(1 - beta), base=2))
  M[M < -threshold] = -threshold
  M[M > threshold] = threshold
  return(M)
}

