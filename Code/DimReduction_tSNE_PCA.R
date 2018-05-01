# This script is for dimensionality reduction and
# examining cell distribution on tSNE and PCA plots.
# It uses the filtered, normalised data in the
# HiPSC_ftGenes_dcvl object

# fontsize setting
fontsize <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 16))
# examine cell distribution on tSNE
# plot
pdf("tSNEplot_overall_distribution_beforeClustering.pdf")
plotTSNE(HiPSC_ftGenes_dcvl, exprs_values = "norm_exprs",
         perplexity = 10) + fontsize
dev.off()
# examine cell distribution on PCA
# plot
pdf("PCAplot_overall_distribution_beforeClustering.pdf")
plotPCA(HiPSC_ftGenes_dcvl, exprs_values = "norm_exprs") +
  fontsize
dev.off()

# run prcomp an write the PCA result
# to @reducedDimension slot
pdf("PCA_dimensional_reduction.pdf")
HiPSC_ftGenes_dcvl_rmMtRb_PCA <- plotPCASCESet(HiPSC_ftGenes_dcvl_rmMtRb,
                                               ntop = 500, ncomponents = 5, exprs_values = "exprs",
                                               return_SCESet = TRUE, scale_features = TRUE,
                                               draw_plot = TRUE, pca_data_input = "exprs",
                                               selected_variables = NULL, theme_size = 10,
                                               legend = "auto")
dev.off()

# quick clustering using Scater
# package, using UMI counts for all
# cells (Not reduced PCA)
HiPSC_ftGenes_dcvl_rmMtRb <- HiPSC_ftGenes_dcvl
clusters <- quickCluster(HiPSC_ftGenes_dcvl_rmMtRb)

# add clusterIDs
HiPSC_ftGenes_dcvl_rmMtRb@phenoData@data$clusterIDs <- clusters

# plot PCA
pdf("plotPCA_quickCluster.pdf")
plotPCA(HiPSC_ftGenes_dcvl_rmMtRb, exprs_values = "norm_exprs",
        colour_by = "clusterIDs") + fontsize
dev.off()

# plot tSNE
pdf("plot_tSNE_quickCluster.pdf")
plotTSNE(HiPSC_ftGenes_dcvl_rmMtRb, exprs_values = "norm_exprs",
         perplexity = 10, colour_by = "clusterIDs") +
  fontsize
dev.off()
