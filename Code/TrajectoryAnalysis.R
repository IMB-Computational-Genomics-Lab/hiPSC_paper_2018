# This script runs two independent analyses: monocle2 and diffusion to compare results with
# the LTS between-cluster transition prediction

library(monocle)
# Read in expression data
load("dat_dcvl_mtrx_T_unLog_minus1_positive_BOC_for_edgeR.Obj")
fulldat <- t
# Read in cluster data
load("my.clusters.Obj")

cluster <- my.clusters

# reformat gene names
fullnames <- rownames(fulldat)
shortnames <- gsub("_.*", "", fullnames)
shortnames <- make.unique(shortnames, sep = "_ADDuniqueSymbol_")
head(fullnames)
head(shortnames)
rownames(fulldat) <- shortnames

# setup phenoDat
cellIDs <- colnames(fulldat)
cluster_cellIDs <- cbind(cluster, cellIDs)
phenoDat <- as.data.frame(cluster_cellIDs)
# row names of the phenoData object should match the column names of the
# expression matrix.
rownames(phenoDat) <- phenoDat$cellIDs

# setup featureDat row names of the featureData object should match row names of
# the expression matrix.
featureDat <- as.data.frame(fullnames)
rownames(featureDat) <- shortnames
colnames(featureDat) <- "FullNames"
# one of the columns of the featureData should be named 'gene_short_name'.
featureDat$gene_short_name <- shortnames

pd <- new("AnnotatedDataFrame", data = phenoDat)
fd <- new("AnnotatedDataFrame", data = featureDat)

head(phenoDat)
head(featureDat)
head(rownames(fulldat))
CarDif <- newCellDataSet(Matrix(fulldat, sparse = T), phenoData = pd, featureData = fd)
dat <- pData(CarDif)

####################################### Prepare DE genes

# combine DE genes for all clusters
DE1vs234 <- read.table("../New_LASSO/result_table_DESeq_1_vs_cluster234_significant10x15K.txt",
                       header = T)
DE1vs234_genes <- gsub("_.*", "", DE1vs234$id)

DE2vs134 <- read.table("../New_LASSO/result_table_DESeq_2_vs_cluster134_significant10x15K.txt",
                       header = T)
DE2vs134_genes <- gsub("_.*", "", DE2vs134$id)

DE3vs124 <- read.table("../New_LASSO/result_table_DESeq_3_vs_cluster234_significant10x15K.txt",
                       header = T)
DE3vs124_genes <- gsub("_.*", "", DE3vs124$id)

DE4vs123 <- read.table("../New_LASSO/result_table_DESeq_4_vs_cluster123_significant10x15K.txt",
                       header = T)
DE4vs123_genes <- gsub("_.*", "", DE4vs123$id)

head(DE1vs234_genes)
head(DE2vs134_genes)
head(DE3vs124_genes)
head(DE4vs123_genes)
Unique_DEgenes <- c(DE1vs234_genes, DE2vs134_genes, DE3vs124_genes, DE4vs123_genes)

DEgenes_idx <- which(row.names(CarDif) %in% Unique_DEgenes)
length(DEgenes_idx)

CarDif_Day_mk <- CarDif[DEgenes_idx, ]

# Set order filter step to specify genes to be used for ordering cells

CarDif_Day_mk <- setOrderingFilter(CarDif_Day_mk, Unique_DEgenes)

CarDif_Day_mk <- estimateSizeFactors(CarDif_Day_mk)

# estimate dispersion (takes a long time to run)

system.time(CarDif_Day_mk <- estimateDispersions(CarDif_Day_mk))

# plot ordered genes
plot_ordering_genes(CarDif_Day_mk)

# dimensionality reduction
CarDif_Day_mk <- reduceDimension(CarDif_Day_mk, max_components = 2, method = "DDRTree")

# order cells
system.time(CarDif_Day_mk <- orderCells(CarDif_Day_mk))

# load cluster information
load("my.clusters.Obj")
cluster <- my.clusters

# plotting monocle trajectory

# color by days
day = 0
pdf(paste0("Monocle_Day", day, "_AllGenes_Trajectory_byClusters.pdf"))
plot_cell_trajectory(dat, color_by = factor(cluster)) + guides(colour = guide_legend(title = "Cluster")) +
  theme(text = element_text(size = 18)) + scale_color_manual(values = c("#8A2022",
  "#CCA47C", "#788E2B", "#2DCBF2"), limits = c("1", "2", "3", "4"))
dev.off()

# color by clusters
cluster2 <- cluster
cluster2[cluster2 != 4] <- 1
plot_cell_trajectory(CarDif_Day_mk, color_by = factor(cluster2))
cluster3 <- cluster
cluster3[cluster3 != 3] <- 1
plot_cell_trajectory(CarDif_Day_mk, color_by = factor(cluster3))

# color by pseudo time
pdf(paste0("Monocle_DEgenes_Trajectory_byPseudotime.pdf"))
plot_cell_trajectory(CarDif_Day_mk, color_by = "Pseudotime") + theme(text = element_text(size = 18))
dev.off()

# color by state
pdf(paste0("Monocle_Day", day, "_KnownMarkers_Trajectory_byState.pdf"))
plot_cell_trajectory(CarDif_Day_mk, color_by = "State")
dev.off()

# color by gene markers
my_genes <- which(row.names(CarDif_Day_mk) %in% c("SMAD3"))
cds_subset <- CarDif_Day_mk[my_genes, ]
pdf(paste0("Monocle_Day", day, "_KnownMarkers_Trajectory_by_Selected_Genes.pdf"))
plot_genes_in_pseudotime(cds_subset, color_by = factor(cluster2)) + guides(colour = guide_legend(title = ""))
dev.off()

# Run diffusion analysis--------------------------------------------------------
library(destiny)

# load normalised count data
diffusionObject <- readRDS("DiffusionMap_allCells_notCPM_notTranspose.RDS")

# run diffusion (better to run on hpc)
dptObject <- DPT(diffusionObject, branching = TRUE)

pdf("DiffusionMap_20PCs_withDPTwithBranching_and_Path.pdf")
plot.DPT(dptObject, divide = 1:3, paths_to = 1:3)
dev.off()

dptObject <- readRDS("Results_dcvl/DiffusionMap_20PCs_notTranspose.RDS")

# load cluster information
load("my.clusters.Obj")
cluster <- my.clusters

# plot with cluster labels
pdf("DiffusionMap_20PCs_colorbyCluster.pdf")
plot(dptObject, col_by = cluster) + guides(colour = guide_legend(title = "Cluster")) +
  theme(text = element_text(size = 18)) + scale_color_manual(values = c("#8A2022",
  "#CCA47C", "#788E2B", "#2DCBF2"), limits = c("1", "2", "3", "4"))

# plot with branches and tips
plot.DPT(dptObject, col_by = "branch", divide = 1:16)
plot.DPT(dptObject, col_tip = "purple", dcs = c(1, 2), col_path = "purple", paths_to = 2)

# plot with differentiation paths
pdf("DiffusionMap_20PCs_withDPTwithBranching_and_Path.pdf")
plot.DPT(dptObject, divide = 1:3, paths_to = 1:3)
