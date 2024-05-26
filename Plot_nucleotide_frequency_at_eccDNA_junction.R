library(ggplot2)
library(ggseqlogo)

#set working directory
setwd('Your/directory/')

##################################################

eccDNA_start <- read.csv("sperm_eccDNA_start_junction_positionFrequency.csv",header=T, stringsAsFactors=F)
eccDNA_end <- read.csv("sperm_eccDNA_end_junction_positionFrequency.csv",header=T, stringsAsFactors=F)

rownames(eccDNA_start) <- eccDNA_start[,1]; eccDNA_start[,1] <- NULL # Move first column to rownames
rownames(eccDNA_end) <- eccDNA_end[,1]; eccDNA_end[,1] <- NULL

eccDNA_start[1:3, 1:4]
eccDNA_end[1:3, 1:4]

eccDNA_start_t = t(eccDNA_start)
eccDNA_end_t = t(eccDNA_end)

eccDNA_start_t[1:3, 1:4]
eccDNA_end_t[1:3, 1:4]

ggseqlogo(eccDNA_start_t, method = 'bits')
ggseqlogo(eccDNA_end_t, method = 'bits')
