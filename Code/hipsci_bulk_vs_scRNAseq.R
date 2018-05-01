#===============================================================================#
# Author: Sam Lukowski
# Date started: 24/1/2017
# Date last updated: 24/1/2017
# 
# About: R script to convert HipSci bulk RNA-seq expression data (ERZ267) to
#  transcript counts and plot against scRNA-seq data. Uses initial processing 
#  steps from tximport vignette. tx2gene.txt supplied.
# 
#===============================================================================#

#===============================================================================#
# COMMANDLINE SECTION (not R) TO OBTAIN HIPSCI DATA
wget -r -np ftp://ftp.sra.ebi.ac.uk/vol1/ERZ267/

#===============================================================================#

# R begins here
source('https://bioconductor.org/biocLite.R')
biocLite('tximport')
install.packages('readr')

library(tximport)
library(readr)

setwd('~/Downloads/hipsci')
my.files <- list.files(recursive=T, pattern = '.tsv')

tx2gene <- read.table('~/Downloads/hipsci/tx2gene.txt', header = T)

txi <- tximport(my.files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
counts <- as.data.frame(txi$counts)
colnames(counts) <- c('HPSI0514i-wiii_2', 'HPSI0314i-cuhk_2', 'HPSI1213i-xuja_2', 'HPSI1113i-wetu_2', 'HPSI1013i-yemz_1', 'HPSI1013i-wuye_2', 'HPSI1013i-jufd_2', 'HPSI0214i-kucg_2', 'HPSI0913i-lise_1', 'HPSI0913i-diku_1', 'HPSI1013i-pamv_3', 'HPSI0114i-fikt_3', 'HPSI0913i-lise_3', 'HPSI0114i-joxm_1', 'HPSI1113i-dons_1', 'HPSI1113i-podx_1', 'HPSI0314i-qonc_1', 'HPSI1113i-hayt_1', 'HPSI0814i-bokz_5', 'HPSI1113i-eofe_1', 'HPSI1113i-qolg_1', 'HPSI0114i-bezi_3', 'HPSI0114i-vabj_3', 'HPSI1013i-sebz_1', 'HPSI0814i-doao_1', 'HPSI0314i-bipt_1', 'HPSI1213i-nusw_2', 'HPSI1013i-jogf_2', 'HPSI0314i-qonc_2', 'HPSI1013i-hiaf_2', 'HPSI0114i-lexy_1', 'HPSI0114i-zapk_3', 'HPSI1213i-pahc_5', 'HPSI1013i-kuxp_1', 'HPSI1213i-babk_2', 'HPSI1113i-qorq_2', 'HPSI0314i-fafq_1', 'HPSI1213i-nekd_1', 'HPSI1013i-wopl_1', 'HPSI0414i-mita_2', 'HPSI1113i-qorq_1', 'HPSI0114i-oevr_3', 'HPSI1113i-bima_2', 'HPSI1013i-pamv_1', 'HPSI0314i-bubh_1', 'HPSI0214i-datg_2', 'HPSI0114i-zoxy_3', 'HPSI0214i-feec_3', 'HPSI0314i-sojd_3', 'HPSI0214i-kehc_2', 'HPSI1113i-bima_1', 'HPSI1113i-qolg_3', 'HPSI0913i-oapg_5', 'HPSI1013i-garx_2', 'HPSI1113i-ieki_3', 'HPSI1013i-cups_3', 'HPSI0214i-heja_2', 'HPSI0114i-eipl_1', 'HPSI0114i-iisa_1', 'HPSI1013i-wuye_3', 'HPSI0114i-rozh_5', 'HPSI0314i-sojd_2', 'HPSI0114i-bezi_1', 'HPSI0214i-rayr_1', 'HPSI0814i-doao_2', 'HPSI0314i-cuhk_1', 'HPSI1113i-ieki_2', 'HPSI1013i-jufd_3', 'HPSI0814i-bokz_6', 'HPSI0214i-eiwy_1', 'HPSI0214i-heth_1')

write.table(counts, 'txiCounts.txt', qu = F, row = T, sep = '\t')


#zcounts <- scale(t(counts), center = T)
#write.table(zcounts, 'zCounts.txt', qu = F, row = T, sep = '\t')


#===============================================================================#
# Plot single cell lasso genes against HipSci bulk RNA-seq



library(dplyr)
library(reshape2)
library(ggplot2)

cluster1 <- c('YBX1', 'GNB2L1', 'NGFRAP1', 'HSPE1', 'MCUR1', 'FTH1', 'PRELID1', 'HMGA1', 'SERF2')
cluster2 <- c('DYNC1LI1', 'SRSF7', 'GSDMD', 'RNMTL1', 'MAD2L1', 'FKBP1A', 'FNTA', 'NUDCD1', 'SAPCD2')
cluster3 <- c('ZDHHC5', 'TENM3', 'NLGN4X', 'COL4A2', 'CASC4', 'CABLES1', 'SESN2', 'PLEC', 'IRS2', 'WRNIP1')
cluster4 <- c('AURKAIP1', 'NDUFS6', 'SIVA1', 'NMT1', 'MAD2L2', 'HSF2', 'SEC11A', 'SLC25A1', 'PTP4A2', 'EZR', 'PGP', 'NRBP1', 'FAM104B', 'ACLY')

hipsci <- read.table('~/Downloads/hipsci/txiCounts.txt', header = T, sep = '\t')
hipsci$geneID <- rownames(hipsci)

hipsci1 <- as.data.frame(hipsci[which(rownames(hipsci) %in% cluster1), ])
hipsci2 <- as.data.frame(hipsci[which(rownames(hipsci) %in% cluster2), ])
hipsci3 <- as.data.frame(hipsci[which(rownames(hipsci) %in% cluster3), ])
hipsci4 <- as.data.frame(hipsci[which(rownames(hipsci) %in% cluster4), ])

hipsci1$type <- 'hipsci'
hipsci2$type <- 'hipsci'
hipsci3$type <- 'hipsci'
hipsci4$type <- 'hipsci'

hipsci1.melt <- melt(hipsci1, varids = 'geneID')
hipsci2.melt <- melt(hipsci2, varids = 'geneID')
hipsci3.melt <- melt(hipsci3, varids = 'geneID')
hipsci4.melt <- melt(hipsci4, varids = 'geneID')

hipsci1.melt$value <- log2(hipsci1.melt$value + 1)
hipsci2.melt$value <- log2(hipsci2.melt$value + 1)
hipsci3.melt$value <- log2(hipsci3.melt$value + 1)
hipsci4.melt$value <- log2(hipsci4.melt$value + 1)



# Let's do the same for the single-cell data!
sc <- readRDS('~/Downloads/Expression_CPM_unLog_minus1_positive.Obj')
sc1 <- as.data.frame(sc[which(rownames(sc) %in% cluster1), ])
sc2 <- as.data.frame(sc[which(rownames(sc) %in% cluster2), ])
sc3 <- as.data.frame(sc[which(rownames(sc) %in% cluster3), ])
sc4 <- as.data.frame(sc[which(rownames(sc) %in% cluster4), ])

sc1$geneID <- rownames(sc1)
sc2$geneID <- rownames(sc2)
sc3$geneID <- rownames(sc3)
sc4$geneID <- rownames(sc4)

sc1$type <- 'singlecell'
sc2$type <- 'singlecell'
sc3$type <- 'singlecell'
sc4$type <- 'singlecell'

sc1.melt <- melt(sc1, varids = 'geneID')
sc2.melt <- melt(sc2, varids = 'geneID')
sc3.melt <- melt(sc3, varids = 'geneID')
sc4.melt <- melt(sc4, varids = 'geneID')

sc1.melt$value <- log2(sc1.melt$value + 1)
sc2.melt$value <- log2(sc2.melt$value + 1)
sc3.melt$value <- log2(sc3.melt$value + 1)
sc4.melt$value <- log2(sc4.melt$value + 1)


# create a combined dataset
combined1 <- rbind(hipsci1.melt, sc1.melt)
combined2 <- rbind(hipsci2.melt, sc2.melt)
combined3 <- rbind(hipsci3.melt, sc3.melt)
combined4 <- rbind(hipsci4.melt, sc4.melt)

# calculate correlations
# correlation between expression trends across technologies:
# remove the geneID and type columns
hipsci1.1 <- hipsci1[, -(72:73)]
hipsci2.1 <- hipsci2[, -(72:73)]
hipsci3.1 <- hipsci3[, -(72:73)]
hipsci4.1 <- hipsci4[, -(72:73)]

sc1.1 <- sc1[, -(18788:18789)]
sc2.1 <- sc2[, -(18788:18789)]
sc3.1 <- sc3[, -(18788:18789)]
sc4.1 <- sc4[, -(18788:18789)]

#calculate row means
hipsci1.mean <- apply(hipsci1.1, 1, mean)
hipsci2.mean <- apply(hipsci2.1, 1, mean)
hipsci3.mean <- apply(hipsci3.1, 1, mean)
hipsci4.mean <- apply(hipsci4.1, 1, mean)

sc1.mean <- apply(sc1.1, 1, mean)
sc2.mean <- apply(sc2.1, 1, mean)
sc3.mean <- apply(sc3.1, 1, mean)
sc4.mean <- apply(sc4.1, 1, mean)


# log2 transform the count data
hipsci1.mean <- log2(hipsci1.mean + 1)
hipsci2.mean <- log2(hipsci2.mean + 1)
hipsci3.mean <- log2(hipsci3.mean + 1)
hipsci4.mean <- log2(hipsci4.mean + 1)

sc1.mean <- log2(sc1.mean + 1)
sc2.mean <- log2(sc2.mean + 1)
sc3.mean <- log2(sc3.mean + 1)
sc4.mean <- log2(sc4.mean + 1)


# sort alphabetically
hipsci1.mean <- hipsci1.mean[order(names(hipsci1.mean))]
hipsci2.mean <- hipsci2.mean[order(names(hipsci2.mean))]
hipsci3.mean <- hipsci3.mean[order(names(hipsci3.mean))]
hipsci4.mean <- hipsci4.mean[order(names(hipsci4.mean))]

sc1.mean <- sc1.mean[order(names(sc1.mean))]
sc2.mean <- sc2.mean[order(names(sc2.mean))]
sc3.mean <- sc3.mean[order(names(sc3.mean))]
sc4.mean <- sc4.mean[order(names(sc4.mean))]


# correlation
cor(hipsci1.mean, sc1.mean) #[1] 0.865902
cor(hipsci2.mean, sc2.mean) 
cor(hipsci3.mean, sc3.mean) 
cor(hipsci4.mean, sc4.mean) 




# generate plots

cols <- c("#E69F00", "#56B4E9") #(orange and blue)

p <- ggplot(combined1, aes(geneID, value, fill = type))
p1 <- p + geom_violin(scale = "width", trim = TRUE, aes(fill = type, geneID, value)) +
	scale_fill_manual(values = cols) +
	xlab('') + 
	ylab('log2 counts') +
	ggtitle('Cluster 1 vs 234') +
	annotate('text', x = 0.8, y = 5, label = paste0('r = ', round(cor(hipsci1.mean, sc1.mean), 3)), size = 3, angle = 90) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(legend.position = 'none') +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
png('~/Desktop/lasso1.png', w = 1500, h = 300, res = 150); print(p1); dev.off()


p <- ggplot(combined2, aes(geneID, value, fill = type))
p2 <- p + geom_violin(scale = "width", trim = TRUE, aes(fill = type, geneID, value)) +
	scale_fill_manual(values = cols) +
	xlab('') + 
	ylab('log2 counts') +
	ggtitle('Cluster 2 vs 134') +
	annotate('text', x = 0.8, y = 4, label = paste0('r = ', round(cor(hipsci2.mean, sc2.mean), 3)), size = 3, angle = 90) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(legend.position = 'none') +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
png('~/Desktop/lasso2.png', w = 1500, h = 300, res = 150); print(p2); dev.off()


p <- ggplot(combined3, aes(geneID, value, fill = type))
p3 <- p + geom_violin(scale = "width", trim = TRUE, aes(fill = type, geneID, value)) +
	scale_fill_manual(values = cols) +
	xlab('') + 
	ylab('log2 counts') +
	ggtitle('Cluster 3 vs 124') +
	annotate('text', x = 0.8, y = 4, label = paste0('r = ', round(cor(hipsci3.mean, sc3.mean), 3)), size = 3, angle = 90, angle = 90) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(legend.position = 'none') +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
png('~/Desktop/lasso3.png', w = 1500, h = 300, res = 150); print(p3); dev.off()


p <- ggplot(combined4, aes(geneID, value, fill = type))
p4 <- p + geom_violin(scale = "width", trim = TRUE, aes(fill = type, geneID, value)) +
	scale_fill_manual(values = cols) +
	xlab('') + 
	ylab('log2 counts') +
	ggtitle('Cluster 4 vs 123') +
	annotate('text', x = 0.8, y = 5, label = paste0('r = ', round(cor(hipsci4.mean, sc4.mean), 3)), size = 3, angle = 90) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(legend.position = 'none') +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	#theme(legend.position = 'bottom', legend.text = element_text(size = 10))
	
png('~/Desktop/lasso4.png', w = 1500, h = 300, res = 150); print(p4); dev.off()


source('~/Documents/working_scripts/R/functions/multiplot.R', chdir = TRUE)

multiplot(p1, p2, p3, p4, cols = 1)
#png('~/Desktop/lasso1-4.png', w = 500, h = 1600, res = 150); multiplot(p1, p2, p3, p4, cols = 1); dev.off()


