#!/usr/bin/env Rscript
# File name: computeCorrelation.R
# Author: Mamie Wang
# Date created: 02/27/2018
# Date modified: 02/28/2018

library(hashmap)

#' Construct a vector of correlation coefficients for each pair of methylation and RNAseq.
#' 
#' @param hash.map a hash map from methylation probe id to gene probe id
#' @param methylation methylation data matrix with row names as methylation probe id
#' @param rnaseq RNAseq data matrix with row names as RNAseq probe id and 
#' same number of columns of the methylation data matrix
#' @return correlation coefficients all possible matching pair of methylation and RNAseq

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
  m.geneid <- cbind(as.character(seq(sum(filtered))), geneid[filtered])
  rownames(m.geneid) <- NULL
  return(list(m.geneid=m.geneid, rnaseq=rnaseq.filtered))
}

#' Preprocess the methylation data for input of MethylationGEPCorrelation
#' 
#' @param methylation Methylation data frame (first column probe id and second 
#' column gene id) 
#' @return a list with probe id - gene id data matrix and data matrix with only
#' methylation level
PreprocessMethylation <- function(methylation) {
  probeid <- methylation$ID
  geneid <- methylation$Gene.Symbol
  notna <- !is.na(geneid)
  methylation.filtered <- methylation[notna, -seq(4)]
  colnames(methylation.filtered) <- sapply(colnames(methylation.filtered),
                                           function(x) tolower(substr(x, 1, 12)))
  rownames(methylation.filtered) <- geneid[notna]
  g.geneid <- cbind(probeid[notna], geneid[notna])
  rownames(g.geneid) <- NULL
  return(list(g.geneid=g.geneid, methylation=methylation.filtered))
}

#' Match the columns between RNAseq and Methylation array
#' 
#' @param rnaseq RNAseq data matrix
#' @param methylation Methylation data matrix
#' @return a list containing RNAseq and methylation with matching columns
matchColumns <- function(rnaseq, methylation) {
  intersection <- intersect(colnames(rnaseq), colnames(methylation))
  intersection.rnaseq <- sapply(intersection, function(x) which(x %in% colnames(rnaseq)))
  intersection.methylation <- sapply(intersection, function(x) which(x %in% colnames(methylation)))
  list(rnaseq=rnaseq[,intersection.rnaseq], 
       methylation=methylation[,intersection.methylation])
}

#' Main function to run the correlation coefficient function
#'
#' @param rnaseq RNAseq data frame
#' @param methylation methylation data frame
#' @return a vector of correlation coefficient between matching pairs of rnaseq 
#' and methylation probes
ComputeCorrelation <- function(rnaseq,  methylation) {
  rnaseq.processed <- PreprocessRNAseq(rnaseq)
  methylation.processed <- PreprocessMethylation(methylation)
  hash.map <- ConstructM2GHashMap(methylation.processed$m.geneid, 
                      rnaseq.processed$g.geneid)
  processed <- matchColumns(rnaseq.processed$rnaseq, 
                            methylation.processed$methylation)
  MethylationGEPCorrelation(hash.map, processed$methylation, processed$rnaseq)
}