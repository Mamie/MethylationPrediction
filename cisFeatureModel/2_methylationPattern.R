# Characterize the methylation pattern
library(ggplot2)
library(FDb.InfiniumMethylation.hg18
load('Beta190.RData')
load('Class190.RData')
rownames(Beta190) <- as.character(Beta190$ID)

# summary statistics of methylation status
Class190 <- data.frame(Class190)
Class190[,seq(2,191)] <- apply(Class190[,seq(2,191)], 2, function(x) as.numeric(as.character(x)))
beta.gt50 <- apply(Class190[,-1], 1, sum)

# histogram of the number of patients with beta > 0.5 for each site
scienceTheme <- theme(
             axis.line = element_line(size=.7, color = "black"),
             legend.position = c(.85,.7),
             text = element_text(size=14))
svg(file='histBeta.svg', width=6, height=4)
ggplot(data=data.frame(x=beta.gt50/190*100)) +
    geom_histogram(aes(x=x)) +
    labs(x=expression('Percent of methylation level > 0.5 (%)'), y='Count') +
    theme_classic(base_size=12, base_family='Helvetica') + 
    scienceTheme
dev.off()
rm(beta.gt50)

# mean beta distribution
beta.mean <- apply(Beta190[,-1], 1, mean)
svg(file='histMeanBeta.svg', width=6, height=4)
ggplot(data=data.frame(x=beta.mean)) +
    geom_histogram(aes(x=x)) +
    labs(x='Mean methylation level', y='Count') +
    theme_classic(base_size=12, base_family='Helvetica') + 
    scienceTheme
dev.off()

mean(beta.mean > 0.7) # 20.89 %
mean(beta.mean < 0.3) # 65.82 %
mean(as.matrix(Class190[,-1])) # 27.23 %
rm(beta.mean)

# Correlation of beta values between neighboring sites
# get genomic coordinate of the probes (hg 27)
annotations <- read.csv('HumanMethylation450_15017482_v1-2.csv', sep=',', stringsAsFactors=F, skip=7, header=T)
row.names(annotations) <- annotations$IlmnID
hm27.probes <- as.character(Class190[,1])
probe.genomicCoordinate <- annotations[hm27.probes, c('CHR', 'MAPINFO')]
probe.genomicCoordinate <- GRanges(seqnames=Rle(paste0('chr', probe.genomicCoordinate$CHR)),
    ranges=IRanges(start=probe.genomicCoordinate$MAPINFO, end=probe.genomicCoordinate$MAPINFO + 1),
    strand=rep('*', dim(probe.genomicCoordinate)[1]))
names(probe.genomicCoordinate) <- hm27.probes
probe.neighbor <- precede(probe.genomicCoordinate)
probe.neighbor <- cbind(hm27.probes, hm27.probes[probe.neighbor])
probe.neighbor <- na.omit(probe.neighbor)
correlation.neighbor <- apply(probe.neighbor, 1, function(x)
                            cor(as.numeric(Beta190[x[1],-1]), as.numeric(Beta190[x[2], -1])))
range(correlation.neighbor, na.rm=T) # -0.9064963  0.9992194

distance.neighbor <- apply(probe.neighbor, 1, function(x)
                            distance(probe.genomicCoordinate[x[1]], probe.genomicCoordinate[x[2]]))
range(distance.neighbor) #  0 28160985
d <- data.frame(x=distance.neighbor, y=abs(correlation.neighbor))
d.spline <- smooth.spline(d)

# plot of correlation of CGI sites methylation level with nearest genomic neighbor on assay with respect to distance
svg(file='neighborCorrelation.svg', width=6, height=4)
ggplot(data=d) +
    geom_point(aes(x=x, y=y), col='slategray1', alpha=0.1) +
    geom_line(data=data.frame(x=d.spline$x, y=d.spline$y), aes(x=x, y=y)) +
    labs(x='CpG site distance (bp)', y='Absolute correlation') +
    scale_x_continuous(limits = c(0, 6000)) +
    theme_classic(base_size=12, base_family='Helvetica') + 
    scienceTheme
dev.off()


save(correlation.neighbor, distance.neighbor, file='correlationsNeighbor.RData')
