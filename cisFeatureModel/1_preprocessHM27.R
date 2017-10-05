library(dplyr)
library(ggplot2)

# Process LAML methylation data
data.HM27 <- read.csv('/rhea/Data/MethylationPrediction/LAML.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt', 
                       sep='\t', header=T, stringsAsFactors=F)
num.columns <- dim(data.HM27)[2]
redundant.idx <- c(seq(7, num.columns, by=4), seq(8, num.columns, by=4), 
                   seq(9, num.columns, by=4))
data.HM27 <- data.HM27[-1,-redundant.idx]
data.HM27[,seq(5)] <- data.HM27[,c(1,3,4,5,2)]
columnnames <- colnames(data.HM27)
columnnames[seq(4)] <- c('ID', 'Gene.Symbol', 'Chromosome', 'Genomic.Coordinate')
colnames(data.HM27) <- columnnames
dim(data.HM27) # 27578 x 198
rm(columnnames, num.columns, redundant.idx)
save(data.HM27, file='AMLHM27.RData')

# Remove rows with chr X/Y/NA, and missing entries
data.HM27 <- data.HM27 %>%
      filter(!(Chromosome %in% c("X", "Y", NA)))
missing.rm <- apply(data.HM27[seq(5, 198)], 1, function(x) sum(is.na(x)) > 0)
data.HM27 <- data.HM27[!missing.rm,]
rm(missing.rm)

# merge with manifest file
data.annot <- read.csv('/rhea/Data/MethylationPrediction/GPL8490-65.txt', sep='\t', skip=38, header=T, stringsAsFactors=F)
data.annot <- data.annot[,-2]
data.HM27 <- merge(data.HM27, data.annot, by='ID', all.x=T)
rm(data.annot)

# Patient gender and age
data.clin <- read.csv('/rhea/Data/MethylationPrediction/LAML.clin.merged.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
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
sum(patient.clin[,1] == patient.id) == 194 # sanity check
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

# histogram plot for 
scienceTheme <- theme(
             axis.line = element_line(size=.7, color = "black"),
             legend.position = c(.85,.7),
             text = element_text(size=14))
beta.hist <- ggplot(data=data.frame(x=c(matrixt))) +
       geom_histogram(aes(x=x)) +
       labs(x="Methylation level", y="Count") +
       theme_classic(base_size = 12, base_family = "Helvetica") +
       scienceTheme
betaAdjusted.hist <- ggplot(data=data.frame(x=c(matrix.control.scale))) +
    geom_histogram(aes(x=x)) +
    labs(x="Normalized methylation level", y="Count") +
    theme_classic(base_size=12, base_family="Helvetica") +
    scienceTheme

svg(file="betaHistPlot.svg", width=6, height=4)
beta.hist
dev.off()
svg(file="normalizedBetaHist.svg", width=6, height=4)
betaAdjusted.hist
dev.off()

# PCA 
project.t <- t(matrix.control.scale)
cor <- cor(project.t)
eigen <- eigen(cor)
loadings <- eigen$vectors
PCA.plot <- ggplot(data=data.frame(PC1=loadings[,1], PC2=loadings[,2])) +
    geom_point(aes(x=PC1, y=PC2)) +
    labs(x='Principal component 1', y='Principal component 2') +
    theme_classic(base_size=12, base_family='Helvetica') +
    scienceTheme
svg(file='PCA.svg', width=6, height=4)
PCA.plot
dev.off()

# remove outliers
outliers <- which(loadings[,1] > -0.065)
patient.clin[outliers,]
project.t <- project.t[,-outliers]
colnames(project.t) <- patient.clin[-outliers, 1]
Beta190 <- data.frame(data.HM27$ID, project.t)
colnames(Beta190)[1] <- 'ID'
rownames(Beta190) <- Beta190[,1]
save(Beta190, file="Beta190.RData")

Class190 <- as.character(Beta190$ID)
for (i in 2:191) {
      class = rep(0, nrow(Beta190))
  class[which(Beta190[,i] >= 0.5)] = 1
    Class190 = cbind(Class190, class)
    colnames(Class190)[i] = paste("V", i - 1, sep = "")
}
colnames(Class190)[1]=c("site")
Class190 <- data.frame(Class190)
rownames(Class190) <- as.character(Beta190[,1])
colnames(Class190) <- colnames(Beta190)
save(Class190, file="Class190.RData")

# Input file links:
# methylation: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LAML/20160128/gdac.broadinstitute.org_LAML.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz
# Illumina HM27 manifest adf file: https://cancergenome.nih.gov/abouttcga/aboutdata/platformdesign
# Clinical data: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LAML/20160128/gdac.broadinstitute.org_LAML.Merge_Clinical.Level_1.2016012800.0.0.tar.gz 
