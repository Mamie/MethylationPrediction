# Helper functions for extract clustering information

library(dplyr)
library(biomaRt)
library(fossil)
library(ggplot2)
library(grid)
library(pheatmap)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")


#' Extract the probe id for given cluster and find corresponding gene name
#'
#' @param cluster the clustering membership for each probe
#' @param k the cluster number
#' @return A data frame containing the probe id and gene name in the given cluster
ExtractClusterProbes <- function(cluster, geneid, k) {
  probe.names <- names(cluster)[cluster==k]
  probe.genes <- geneid %>%
    dplyr::filter(ID %in% probe.names) %>%
    dplyr::select(ID, Gene.Symbol)
  # split fields with multiple gene symbol as separate rows
  n.rows <- dim(probe.genes)[1]
  ids <- c()
  genes <- c()
  for (i in seq(n.rows)) {
      id <- probe.genes[i,  1]
      gene <- probe.genes[i, 2]
      gene <- strsplit(gene, ',')[[1]]
      n.gene <- length(gene)
      if (n.gene > 1) {
          ids <- c(ids, rep(id, n.gene))
          genes <- c(genes, gene)
      } else {
          ids <- c(ids, id)
          genes <- c(genes, gene)
      }
  }
  print(paste('The number of probes for cluster', k, 'is', n.rows, 'and the number of genes is', length(ids)))
  return(data.frame(id=ids, gene=genes))
}

#' Merge the hgnc names with the hgnc description, NA if not found
#'
#' @param geneid a data frame containing the gene and probe id
#' @return a data frame with gene description from hgnc if available
MergeHGNCDescription <- function(geneid) {
    genes.info <- getBM(c("hgnc_symbol", "description"), "hgnc_symbol", geneid$gene, mart=ensembl)
    if (dim(genes.info)[1] == 0) return(NULL)
    geneid %>%
      left_join(genes.info, by=c('gene'='hgnc_symbol'))
}

#' Write csv file for data frame
#'
#' @param data the data frame to be written
#' @param path the writing path
#' @param None (side effect: write data frame to csv)
WriteCSV <- function(data, path) {
    write.table(data, file=path, sep=',', quote=F, 
               row.names=F, col.names=T)
    print(paste0('Wrote csv to ', path))
}

#' Wrapper function to extract cluster probes and merge description and write csv
#'
#' @param cluster a vector indicating cluster membership
#' @param geneid a data frame containing mapping information from probe id to hgnc gene symbol
#' @param folderpath the folder to write to write the result csv
WriteClusterInfo <- function(cluster, geneid, folderpath) {
    clusters <- unique(cluster)
    for (k in clusters) {
        clusterk <- ExtractClusterProbes(cluster, geneid, k)
        mergedk <- MergeHGNCDescription(clusterk)
        path <- paste0(folderpath, '/cluster', k, '.csv')
        WriteCSV(mergedk, path)
    }
}

#' Wrapper function to extract the predictor information of the clusters and merge with description and write csv
#'
#' @param folder the path to the folder containing RData of the predictors
#' @param k the number of clusters
#' @return None (side effect: write the information to the same folder)
WritePredictorInfo <- function(folder, k) {
    for (i in seq(k)) {
        base.name <- paste0('/cluster', i, '.RData')
        path <- paste0(folder, base.name)
        out.path <- paste0(folder, '/cluster', i, 'Predictor.csv')
        if(file.exists(path)) {
            load(path)
            if(dim(coef)[1] == 0) {
                print(paste("No predictor selected for cluster", i))
                next
            }
            else {
                merged <- MergeHGNCDescription(coef)
                if (!is.null(merged)) WriteCSV(merged, out.path)
            }
        } else {
            next
        }
    }
}

  
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
# example usage: save_pheatmap_pdf(xx, "test.pdf") where ww is the pheatmap object