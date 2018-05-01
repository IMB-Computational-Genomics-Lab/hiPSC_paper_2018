# This script is to preprocess data from 10X cellRanger
# mapped outputs (merged samples, already normalised
# between samples by cellRanger subsampling procedure).
# 1) It reports the quality, and filters cells and genes outliers
# 2) It performs cell to cell normalisation
#This script can be run automated on HPC

argv <- commandArgs(TRUE)

library(cellrangerRkit)
library(dplyr)
library(Matrix)
library(gdata)
library(R.utils)
library(scran)
library(limSolve)
library(reshape2)
library(dynamicTreeCut)
library(scater)
library(edgeR)

path = argv[1]
# path to parent folder containing
# cellRanger output, e.g. SI-3A-A3....

setwd(path)
print(c("PBS path:", path, "\n"))
dir.create("SC_Cellranger_BOC_5samples")

path_analysis = paste0(path, "SC_Cellranger_BOC_5samples")
print(c("Output path:", path_analysis,
        "\n"))
setwd(path_analysis)
# load expression data
# (exprs_merged_JMpl object)
load(file = paste0(path, "exprs_merged_day0HiSPC_raw.Obj"))

# create a scater object
sce_HiPSC <- newSCESet(exprsData = exprs_merged_JMpl,
                       countData = exprs_merged_JMpl)

# Find mitochondrial genes in our
# HiPSC dataset
is.mito_HiPSC <- grepl("^MT-", rownames(sce_HiPSC))

# Find ribosomal genes
is.ribosomal_HiPSC <- grepl("^RPS|^RPL",
                            rownames(sce_HiPSC))
# Can also add cell cyles here and
# perform similar QC

# Calculate QC matrix for each cell,
# stored in pData os the SCEset
sce_HiPSC <- calculateQCMetrics(sce_HiPSC,
                                feature_controls = list(Rb = is.ribosomal_HiPSC,
                                                        Mt = is.mito_HiPSC))

# Plots total genes and library sizes
png("PlotQC_Mito_Ribo.png", w = 2000,
    h = 2000, res = 400)
par(mfrow = c(2, 2), cex = 1.2)
hist(sce_HiPSC$total_counts/1e+06, xlab = "Library sizes (millions)",
     main = "", breaks = 20, col = "grey80",
     ylab = "Number of cells")
hist(sce_HiPSC$total_features, xlab = "Number of expressed genes",
     main = "", breaks = 20, col = "grey80",
     ylab = "Number of cells")

# plots reads mapped to Mt genes
hist(sce_HiPSC$pct_counts_feature_controls_Mt,
     xlab = "Mitochondrial proportion (%)",
     ylab = "Number of cells", breaks = 20,
     main = "", col = "grey80")

# plots reads mapped to Rb genes
hist(sce_HiPSC$pct_counts_feature_controls_Rb,
     xlab = "Ribosomal proportion (%)",
     ylab = "Number of cells", breaks = 20,
     main = "", col = "grey80")
dev.off()
png("AverageCount.png", w = 2000, h = 2000,
    res = 400)

# examine expression of log-means
# across all genes
ave.counts <- rowMeans(counts(sce_HiPSC))
hist(log10(ave.counts), breaks = 100,
     main = "", col = "grey80", xlab = expression(Log[10] ~
                                                    "average count"))

# Plot number of top genes
fontsize <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 16))
png("TopGenes_gotMapped.png", w = 2000,
    h = 2000, res = 400)
plotQC(sce_HiPSC, type = "highest-expression",
       n = 30) + fontsize
dev.off()

# Plot number of cells
numcells <- nexprs(sce_HiPSC, byrow = TRUE)
png("AverageCount_SmoothScatter.png",
    w = 2000, h = 2000, res = 400)
smoothScatter(log10(ave.counts), numcells,
              xlab = expression(Log[10] ~ "average count"),
              ylab = "Number of expressing cells")
dev.off()

# PCA plot to check for potential cell
# outliers according to PC1 and PC2
# (based on general cell data)
png("PCA_cellOutliers.png", w = 2000,
    h = 2000, res = 400)
plotPCA(sce_HiPSC, pca_data_input = "pdata") +
  fontsize
dev.off()

############## Additional data cleaning
############## steps####################################

# remove cells with low expression or
# low number of genes (lower than 3
# median absolute deviation of
# log(library size))
libsize.drop_HiPSC <- isOutlier(sce_HiPSC$total_counts,
                                nmads = 3, type = "lower", log = TRUE)
feature.drop_HiPSC <- isOutlier(sce_HiPSC$total_features,
                                nmads = 3, type = "lower", log = TRUE)

# remove cells with high percent of
# reads mapped to Mt genes (possibly
# dead cells) (higher than 3 median
# absolute deviation)
mito.drop_HiPSC <- isOutlier(sce_HiPSC$pct_counts_feature_controls_Mt,
                             nmads = 3, type = "higher")
ribo.drop_HiPSC <- isOutlier(sce_HiPSC$pct_counts_feature_controls_Rb,
                             nmads = 3, type = "higher")
sce_HiPSC <- sce_HiPSC[, !(libsize.drop_HiPSC |
                             feature.drop_HiPSC | mito.drop_HiPSC |
                             ribo.drop_HiPSC)]

# Remove cells by Mt gene further
mito.drop0.2 <- sce_HiPSC$pct_counts_feature_controls_Mt <=
  20
mito.remove0.2 <- sce_HiPSC$pct_counts_feature_controls_Mt >
  20
sce_HiPSC <- sce_HiPSC[, mito.drop0.2]

# Remove cells by Rb gene further
ribo.drop0.5 <- sce_HiPSC$pct_counts_feature_controls_Rb <=
  50
ribo.remove0.5 <- sce_HiPSC$pct_counts_feature_controls_Rb >
  50
sce_HiPSC <- sce_HiPSC[, ribo.drop0.5]

# check number of cells that genes
# express
png("Cells_expresing_genes.png", w = 2000,
    h = 2000, res = 400)
numcells <- nexprs(sce_HiPSC, byrow = TRUE)
hist(log2(numcells), xlab = "Log2 number of cells expressing the gene",
     ylab = "Number of genes", main = "Number of cells a gene was detected")
dev.off()

# remove genes expressed in fewer than
# 1% of total cells
numcells <- nexprs(sce_HiPSC, byrow = TRUE)
genes.keep <- numcells >= 22
genes.remove <- numcells < 22
sce_HiPSC_ftGenes <- sce_HiPSC[genes.keep,
                               ]

# Plot number genes after data
# filtering
png("TopGenes_gotMapped_PostDataCleaning.png")
plotQC(sce_HiPSC_ftGenes, type = "highest-expression",
       n = 30) + fontsize
dev.off()

# check how many cells are removed
datRemove <- data.frame(ByLibSize = sum(libsize.drop_HiPSC),
                        ByFeature = sum(feature.drop_HiPSC),
                        GeneRemovedByCell = sum(genes.remove),
                        ByMito1 = sum(mito.drop_HiPSC), ByMito0.2 = sum(mito.remove0.2),
                        byRibo1 = sum(ribo.drop_HiPSC), ByRibo0.5 = sum(ribo.remove0.5),
                        CellRemaining = ncol(sce_HiPSC_ftGenes),
                        GeneRemaining = sum(genes.keep))

# note Mt and Rb genes were NOT at
# this stage

write.table(datRemove, "number_cells_genes_removed.txt",
            quote = F, row.names = F, col.names = T,
            sep = "\t")
save(sce_HiPSC_ftGenes, file = "HiPSC_ftGenes_ReadyFor_ComputeSumFactors_BOC.Obj")

########## Remove Mt genes and Rb genes before
########## clustering########################
mito_ftGenes <- grep("^MT-", rownames(sce_HiPSC_ftGenes))

ribosomal_ftGenes <- grep("^RPL|^RPS",
                          rownames(sce_HiPSC_ftGenes))

mito_ribo <- c(mito_ftGenes, ribosomal_ftGenes)

sce_HiPSC_ftGenes_rmMtRb <- sce_HiPSC_ftGenes[-mito_ribo,
                                              ]

######### Normalization by deconvolution
######### method##################################
library(scran)
library(limSolve)

# Caution: takes a long time on
# laptop, better run it in cluster

sce_HiPSC_ftGenes_cp <- computeSumFactors(sce_HiPSC_ftGenes_rmMtRb,
                                          sizes = c(40, 60, 80, 100), positive = T)

# Remove zero size factors by
# converting them to the minimum size
# factor value:
Zero_sizefactor <- which(sce_HiPSC_ftGenes_cp@phenoData@data$size_factor ==
                           0)
min_size <- min(sce_HiPSC_ftGenes_cp@phenoData@data$size_factor[-Zero_sizefactor])
sce_HiPSC_ftGenes_cp@phenoData@data$size_factor[Zero_sizefactor] <- min_size

# plot the normalised data
png("sizeFactors_normalized.png")
plot(sizeFactors(sce_HiPSC_ftGenes_cp),
     sce_HiPSC_ftGenes_cp$total_counts/1e+06,
     log = "xy", ylab = "Library size (millions)",
     xlab = "Size factor")
dev.off()
# save the object before normalisation
save(sce_HiPSC_ftGenes_cp, file = "HiPSC_ftGenes_ComputeSizeFactor_before_dcvl.Obj")

# perform the KEY NORMALIZATION STEP
sce_HiPSC_ftGenes_dcvl <- normalize(sce_HiPSC_ftGenesi_cp)
# save the object after normalisation
save(sce_HiPSC_ftGenes_dcvl, file = "HiPSC_ftGenes_dcvl.Obj")
