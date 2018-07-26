# Construct feature matrix for user tcga.ab.2803

load('/rhea/Data/MethylationPrediction/feature450K.RData')
load('Beta190.RData')
load('Class190.RData')
library(GenomicRanges)

rownames(shuffle) <- as.character(shuffle$site)
Class190[, seq(2,191)] <- apply(Class190[, seq(2, 191)], 2, function(x) as.character(as.numeric(x)))
annotations <- read.csv('/rhea/Data/MethylationPrediction/HumanMethylation450_15017482_v1-2.csv', sep=',', stringsAsFactors=F, skip=7, header=T)
rownames(annotations) <- annotations$IlmnID


hm27.probes <- as.character(Class190[,1])
features <- data.frame(site=hm27.probes)
rownames(features) <- features$site
features$beta <- Beta190[hm27.probes, 3] # tcga.ab.2803
features$class <- Class190[hm27.probes, 3]

probe.genomicCoordinate <- annotations[hm27.probes, c('CHR', 'MAPINFO')]
probe.genomicCoordinate <- GRanges(seqnames=Rle(paste0('chr', probe.genomicCoordinate$CHR)),
                                       ranges=IRanges(start=probe.genomicCoordinate$MAPINFO, end=probe.genomicCoordinate$MAPINFO + 1),
                                           strand=rep('*', dim(probe.genomicCoordinate)[1]))
names(probe.genomicCoordinate) <- hm27.probes
probe.up <- follow(probe.genomicCoordinate)
probe.up <- hm27.probes[probe.up]
features$site_up <- probe.up

probe.down <- precede(probe.genomicCoordinate)
probe.down <- hm27.probes[probe.down]
features$site_down <- probe.down
features <- na.omit(features)
features$site <- as.character(features$site)

features$beta_up <- Beta190[features$site_up, 3] 
features$class_up <- Class190[features$site_up, 3]
features$beta_down <- Beta190[features$site_down, 3]
features$class_down <- Class190[features$site_down, 3]
features$CHR <- as.character(annotations[features$site, 'CHR'])
features$MAPINFO <- annotations[features$site, 'MAPINFO']
features$MAPINFO_up <- annotations[features$site_up, 'MAPINFO']
features$MAPINFO_down <- annotations[features$site_down, 'MAPINFO']
features$dist_up <- features$MAPINFO - features$MAPINFO_up - 1
features$dist_down <- features$MAPINFO_down - features$MAPINFO - 1
features$PROBE_SNPS <- as.character(as.numeric(annotations[features$site, 'Probe_SNPs'] != ''))
features$PROBE_SNPS_10 <- as.character(as.numeric(annotations[features$site, 'Probe_SNPs_10'] != ''))
# Ref gene group has seven categories
# "1stExon" "5'UTR"   "Body"    "TSS1500" "TSS200"  ""        "3'UTR"
features$island <- as.character(as.numeric(grepl('Island', annotations[features$site, 'Relation_to_UCSC_CpG_Island'])))
features$shore <- as.character(as.numeric(grepl('Shore', annotations[features$site, 'Relation_to_UCSC_CpG_Island'])))
features$shelf <- as.character(as.numeric(grepl('Shelf', annotations[features$site, 'Relation_to_UCSC_CpG_Island'])))
features$nonisland <- as.character(as.numeric(annotations[features$site, 'Relation_to_UCSC_CpG_Island'] == ''))
features$firstExon <- as.character(as.numeric(grepl('1stExon', annotations[features$site, 'UCSC_RefGene_Group'])))
features$fivePrimeUTR <- as.character(as.numeric(grepl("5'UTR", annotations[features$site, 'UCSC_RefGene_Group'])))
features$threePrimeUTR <- as.character(as.numeric(grepl("3'UTR", annotations[features$site, 'UCSC_RefGene_Group'])))
features$geneBody <- as.character(as.numeric(grepl("Body", annotations[features$site, 'UCSC_RefGene_Group'])))
features$TSS1500 <- as.character(as.numeric(grepl("TSS1500", annotations[features$site, 'UCSC_RefGene_Group'])))
features$TSS200 <- as.character(as.numeric(grepl("TSS200", annotations[features$site, 'UCSC_RefGene_Group'])))
features$other <- as.character(as.numeric(annotations[features$site, 'UCSC_RefGene_Group']==''))
geneGroup.site <- as.matrix(features[features$site, seq(22, 27)])
geneGroup.site[,] <- geneGroup.site %in% c('1')
geneGroup.site <- sapply(as.data.frame(geneGroup.site), as.logical)
geneGroup.siteUp <- as.matrix(features[features$site_up, seq(22, 27)])
geneGroup.siteUp[,] <- geneGroup.siteUp %in% c('1')
geneGroup.siteUp <- sapply(as.data.frame(geneGroup.siteUp), as.logical)
geneGroup.siteDown <- as.matrix(features[features$site_down, seq(22, 27)])
geneGroup.siteDown[,] <- geneGroup.siteDown %in% c('1')
geneGroup.siteDown <- sapply(as.data.frame(geneGroup.siteDown), as.logical)
features$gene_group_same_region_up <- as.character(as.numeric(apply(geneGroup.site & geneGroup.siteUp, 1, function(x) sum(x) > 0)))
features$gene_group_same_region_down <- as.character(as.numeric(apply(geneGroup.site & geneGroup.siteDown, 1, function(x) sum(x) > 0)))
island.info <- cbind(annotations[features$site, 'UCSC_CpG_Islands_Name'], annotations[features$site_down, 'UCSC_CpG_Islands_Name'], annotations[features$site_up, 'UCSC_CpG_Islands_Name'])
features$island_same_region_up <- as.character(as.numeric(apply(island.info, 1, function(x) !(any(x[c(1,3)] %in% c(''))) & x[1]==x[3])))
features$island_same_region_down <- as.character(as.numeric(apply(island.info, 1, function(x) !(any(x[c(1,2)] %in% c(''))) & x[1]==x[2])))
features$DHS <- as.character(as.numeric(!is.na(annotations[features$site, 'DHS'])))

intersect.probes <- intersect(features$site, shuffle$site)
features <- cbind(features[intersect.probes,], shuffle[intersect.probes, seq(22, 103)])
features <- cbind(features, shuffle[intersect.probes, c(104, seq(116, 140))])
dim(features) # 17570 x 141
save(features, file='feature2803.RData')
