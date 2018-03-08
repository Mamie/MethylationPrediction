#!/usr/bin/env Rscript
# File name: compareClustering.R
# Author: Mamie Wang
# Date created: 03/08/2018
# Date modified: 03/08/2018

library(fossil)

load('../../data/180307_M/HM450/resultsM/clusterInfo.RData')
HM450.clusterM <- clustering$cluster
ordering <- order(names(clustering$cluster))
HM450.clusterM <- clustering$cluster[ordering]

load('../../data/180307_M/HM450/resultsB/clusterInfo.RData')
HM450.clusterB <- clustering$cluster
ordering <- order(names(clustering$cluster))
HM450.clusterB <- clustering$cluster[ordering]
HM450.rand <- rand.index(HM450.clusterM, HM450.clusterB)
print(paste('Rand index of HM450 clustering is', HM450.rand))

load('../../data/180307_M/HM27/resultsM/clusterInfo.RData')
HM27.clusterM <- clustering$cluster
ordering <- order(names(clustering$cluster))
HM27.clusterM <- clustering$cluster[ordering]

load('../../data/180307_M/HM27/resultsB/clusterInfo.RData')
HM27.clusterB <- clustering$cluster
ordering <- order(names(clustering$cluster))
HM27.clusterB <- clustering$cluster[ordering]
HM27.rand <- rand.index(HM27.clusterM, HM27.clusterB)
print(paste('Rand index of HM27 clustering is', HM27.rand))


