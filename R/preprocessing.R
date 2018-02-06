# File name: preprocessing.R
# Author: Mamie Wang
# Date created: 02/05/2018
# Date last modified: 02/05/2018
# R version: 3.4.2
# Preprocessing RNA-seq: 
#   - remove NA (saved as data/processed/RNAseq.RData)
#   - filter out low variance genes (sd < 0.5)
#   - get annotation for probes and filter for those with terms 'transcrip' or 
#     'methylation'
# Preprocess methylation:
#   - remove NA, chrX or chrY
#   - control for age, gender
#   - filter for low variance probes
#   - hclust use Euclidean distance and Ward method


library(dplyr)
library(tidyr)
library(mygene)
library(ggplot2)


methClusterNums = 20

# Process RNA-seq data: remove NA
print('Reading in mRNA-seq data...')
RNAseq <- read.csv('../data/geneExpression/LAML.uncv2.mRNAseq_RSEM_normalized_log2.txt',
                   sep='\t', header=T, stringsAsFactors=F)
print('Preprocessing mRNA-seq data...')
RNAseq <- na.omit(RNAseq)
FormatPatientID <- function(x) substr(tolower(x), 1, 12)
colnames(RNAseq) <- sapply(colnames(RNAseq), FormatPatientID)
save(RNAseq, file='../data/processed/RNAseq.RData')

# Process methylation data
print('Reading in methylation data...')
data.HM27 <- read.csv('../data/methylationHM27/LAML.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt', 
                      sep='\t', header=T, stringsAsFactors=F)
print('Preprocessing methylation data...')
num.columns <- dim(data.HM27)[2]
redundant.idx <- c(seq(7, num.columns, by=4), seq(8, num.columns, by=4), 
                   seq(9, num.columns, by=4))
data.HM27 <- data.HM27[-1,-redundant.idx]
data.HM27[,seq(5)] <- data.HM27[,c(1,3,4,5,2)]
columnnames <- colnames(data.HM27)
columnnames[seq(4)] <- c('ID', 'Gene.Symbol', 'Chromosome', 'Genomic.Coordinate')
colnames(data.HM27) <- columnnames
rm(columnnames, num.columns, redundant.idx)

# Remove rows with chr X/Y/NA, and missing entries
data.HM27 <- data.HM27 %>%
  filter(!(Chromosome %in% c("X", "Y", NA)))
missing.rm <- apply(data.HM27[seq(5, 198)], 1, function(x) sum(is.na(x)) > 0)
data.HM27 <- data.HM27[!missing.rm,]
rm(missing.rm)

# merge with probe information file
data.annot <- read.csv('../data/GPL8490-65.txt', sep='\t', skip=38, header=T, stringsAsFactors=F)
data.annot <- data.annot[,-2]
data.HM27 <- merge(data.HM27, data.annot, by='ID', all.x=T)
probeGene <- select(data.HM27, ID, Gene.Symbol)
rm(data.annot)

# Patient gender and age
data.clin <- read.csv('../data/clinical/LAML.clin.merged.txt', sep='\t', header=T, 
                      row.names=1, stringsAsFactors=F)
patient.vars <- c('patient.bcr_patient_barcode', 'patient.days_to_birth', 'patient.gender')
data.clin <- data.clin[patient.vars, ]
data.clin <- t(data.clin)
colnames(data.clin) <- c('ID', 'Age', 'Gender')
data.clin[, 'ID'] <- sapply(data.clin[, 'ID'], function(x) gsub('-', '.', x, fixed=T))
rm(patient.vars)

# Normalization for gender and age
matrixt <- as.matrix(t(data.HM27[, seq(5, 198)]))
class(matrixt) <- "double"
patient.id <- sapply(rownames(matrixt), function(x) tolower(substr(x, 1, 12)))
matrixt <- matrixt[order(patient.id), ]
patient.clin <- data.clin[data.clin[,1] %in% patient.id,]
patient.clin <- patient.clin[order(patient.clin[,1]),]
age <- -as.numeric(patient.clin[,'Age'])/365
gender <- as.factor(patient.clin[, 'Gender'])
rm(patient.id)

design.matrix <- model.matrix(~ age + gender)
lm.fit <- lm(matrixt ~ . + 0, data=as.data.frame(design.matrix))
matrix.control <- matrix(nrow = nrow(lm.fit$residuals), ncol = ncol(lm.fit$residuals))
for (i in 1:dim(matrixt)[2]){
  matrix.control[,i] <- c(lm.fit$residuals[,i]+lm.fit$coefficients[1,i])
}

matrix.control.scale <- matrix.control
matrix.control.scale[which(matrix.control > 1)] <- 1
matrix.control.scale[which(matrix.control.scale < 0)] <- 0

site.gt1 <- which(matrix.control> 1)
site.sl0 <- which(matrix.control< 0)

# PCA 
project.t <- t(matrix.control.scale)
cor <- cor(project.t)
eigen <- eigen(cor)
loadings <- eigen$vectors

# visual insepction of the PCA plot and remove outliers
outliers <- which(loadings[,1] > -0.065)
# patient.clin[outliers,]
project.t <- project.t[,-outliers]
colnames(project.t) <- patient.clin[-outliers, 1]
Beta190 <- data.frame(data.HM27$ID, project.t)
colnames(Beta190)[1] <- 'ID'
rownames(Beta190) <- Beta190[,1]
save(Beta190, file="../data/processed/Beta190.RData")
rm(data.HM27, outliers, project.t, cor, eigen, loadings, site.gt1, 
   matrix.control.scale, lm.fit, design.matrix, age, gender, site.sl0,
   FormatPatientID, data.clin, matrix.control, matrixt, patient.clin)

# Get matching patients
print('Intersecting patients between the two data...')
patients.intersect <- intersect(colnames(RNAseq), colnames(Beta190))
M <- Beta190[, colnames(Beta190) %in% patients.intersect]
M <- data.matrix(M)


# Compute the covariance matrix of the methylation probes
print('Computing covariance matrix of methylation level... It will take several minutes')
M.cov <- cov(t(M), method='pearson')
M.absCov <- abs(M.cov)
M.sumAbsCov <- apply(M.absCov, 1, sum)

# kmeans clustering on the sum of absolute covariance. Plot the clustering.
set.seed(1)
cluster.sumAbsCov <- kmeans(M.sumAbsCov, 2)
png(file='../data/figures/filterMethylationKmeans.png')
plot(c(M.sumAbsCov[cluster.sumAbsCov[[1]]==1], M.sumAbsCov[cluster.sumAbsCov[[1]]==2]), 
     ylab='Sum of absolute covariance', xlab='methylation probes')
dev.off()

# Remove the probes in cluster with small sum of covariance
print('Filtering low variance methylation probes...')
M.filtered <- M[cluster.sumAbsCov[[1]]==2, ]
save(M.filtered, file='../data/processed/MethylationFiltered.RData')
rm(M, M.cov, M.absCov, M.sumAbsCov)

# filter out low variance mRNA-seq probes
G.genes <- RNAseq[,1]
G <- RNAseq[, colnames(RNAseq) %in% patients.intersect]
G <- data.matrix(G)
G.std <- apply(G, 1, sd)
rm(RNAseq, patients.intersect)

mask <- G.std > 0.5
G.filtered <- cbind(data.frame(genes=as.character(G.genes[mask])), G[mask,])
save(G.filtered, file='../data/processed/filteredGEP.RData')
rm(mask, G)

# Get annotation for each gene 
print('Getting annotations for mRNA-seq genes')
genes.GEP <- unique(sapply(as.character(G.filtered[,1]), 
                    function(x) strsplit(x, '|', fixed=T)[[1]][1]))
res <- queryMany(genes.GEP, 
                 scopes='symbol', fields=c('go'), species='human',
                 returnall=T)
terms <- c()
for (i in seq(length(res$response$query))) {
  if (is.na(res$response$notfound[i])) {
    query <- res$response$query[i]
    termConcat <- paste(as.character(res$response$go.BP[[i]]$term), collapse=', ')
    terms <- rbind(terms, c(query, termConcat))
  }
}

# Filter for annotation with keywords "methylation", "transcrip". 
print('Filtering gene probes with keywords related to methylation...')
terms.cleaned <- data.frame(name=terms[,1], annotation=terms[,2]) %>%
  dplyr::group_by(name) %>%
  dplyr::mutate(annotations = paste0(annotation, collapse = ", ")) %>%
  dplyr::select(-annotation) %>%
  dplyr::distinct() 
terms.filtered <- terms.cleaned %>%
  dplyr::filter(grepl('methyl|transcrip', annotations)) # 1629 left
save(terms.filtered, file='../data/processed/transcriptionRelatedGenes.RData')
rm(terms.cleaned, res, terms, query, G.std, termConcat, G.genes)

# Hierachical clustering to define methylation clusters
print('Hierarchical clustering for methylation')
distances <- dist(M.filtered) # euclidean distance
hclust.M <- hclust(distances, method='ward.D2')

# select probe with genes that pass the keyword filtration
probeGene$ID <- as.character(probeGene$ID)
probeNames <- data.frame(probe=rownames(M.filtered), cluster=cutree(hclust.M, methClusterNums)) %>% 
  dplyr::left_join(probeGene, by=c('probe'='ID'))

geneSymbols <- strsplit(probeNames$Gene.Symbol, ";")
max_len <- max(sapply(geneSymbols, length))
for(i in 1:max_len){
  probeNames[,paste0("Gene.Symbol", i)] <- sapply(geneSymbols, "[", i)
}

genes <- probeNames %>%
  dplyr::select(-Gene.Symbol) %>%
  tidyr::gather(Gene.Symbol, Name, Gene.Symbol1:Gene.Symbol21) %>%
  dplyr::filter(!is.na(Name)) 

save(genes, hclust.M, file='../data/processed/clusterTCGA.RData')
