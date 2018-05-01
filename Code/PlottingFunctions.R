# This script produces violin and tSNE plots in Figure 2

# tSNE plots--------------------------------------------------------------------
# load tSNE 2 dimensional daaset
load("HiPSC_ftGenes_dcvl_tSNE_2D_scaterObj2.Obj")
tSNE_2D_scaterObj2 <- HiPSC_ftGenes_dcvl_2D@reducedDimension

# load expression dataset
ori_dat_cpm <- readRDS(file = "Expression_CPM_unLog_minus1_positive.Obj")

# load tSNE 3D dataset
dat3d <- readRDS("tSNE_3D_data_with_cluster_batchInfo.Obj")

tSNE_2D_scaterObj2_merged <- cbind(tSNE_2D_scaterObj2, dat3d[, 4:5])
colnames(tSNE_2D_scaterObj2_merged) <- c("tSNE1", "tSNE2", "batches", "cluster")


# plot cells on tSNE with cluster colors
p_cluster <- ggplot(tSNE_2D_scaterObj2_merged, aes(tSNE1, tSNE2, color = factor(cluster))) +
  geom_point(alpha = 0.8) + scale_color_manual(values = c("#E3BA22", "#E6842A",
                                                          "#16687F", "#8E6D8A"), "Cluster")

p_cluster <- p_cluster + theme_bw() + theme(axis.text = element_text(size = 18)) +
  theme() + theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 18)) + guides(color = FALSE)

p_cluster

# prepare tSNE data input (cells with expression values)

tSNEdat = as.data.frame(tSNE_2D_scaterObj2)
colnames(tSNEdat) <- c("tSNE1", "tSNE2")
tSNEdat$cluster <- tSNE_2D_scaterObj2_merged$cluster
tSNEdat$size <- 0.1

AllGeneNames = row.names(ori_dat_cpm)
AllGeneNames <- gsub("_.*", "", AllGeneNames)
SelectedGeneName = "OTX2"
ExpressionMat = ori_dat_cpm

tSNE_gene <- function(AllGeneNames, SelectedGeneName, ExpressionMat, tSNEdat) {
  expression = ExpressionMat
  gene_idx <- which(AllGeneNames == SelectedGeneName)
  GeneExpression <- expression[gene_idx, ]
  positive_index <- which(GeneExpression > 0)
  tSNEdat$size[positive_index] <- 0.4
  Log2GeneExpression <- log2(GeneExpression + 2)
  p <- ggplot(tSNEdat, aes(tSNE1, tSNE2)) + theme_bw() + geom_point(aes(size = size,
                                                                        color = Log2GeneExpression), alpha = 0.4) + scale_colour_gradient2(low = "white",
                                                                                                                                           mid = "#208eb7", high = "#ec4b18") + guides(color = FALSE) + guides(size = FALSE)
  p <- p + theme(text = element_text(size = 12)) + theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank())
  p
}

genes <- c("SOX2", "KLF4", "OTX2", "POU5F1", "LEFTY2", "NANOG", "DNMT3B", "DPPA5",
           "DPPA2", "NODAL", "UTF1")

pathsave <- "./"

tSNE_gene(AllGeneNames = AllGeneNames, SelectedGeneName = "NANOG", ExpressionMat = ori_dat_cpm,
          tSNEdat = tSNEdat)

ggsave(paste0(pathsave, "NANOG_tSNE_Expression.png"), height = 78, width = 60, unit = "mm",
       device = "png", dpi = 300)

# Violin plot--------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(reshape2)

# read exprs converted to CPM
ori_dat_cpm <- readRDS(file = "Expression_CPM_unLog_minus1_positive.Obj")
load("my.clusters.Obj")

# plot log2 CPM values
Violin_single = function(gene, data, cell.ident, ylab.max = 12, size.x.use = 8, size.y.use = 8,
                         size.title.use = 12, adjust.use = 1, size.use = 0.1) {
  data.use = data[which(rownames(data) %in% gene), ]
  data.use <- as.data.frame(data.use)
  data.use$gene <- row.names(data.use)
  data.melt = melt(data.use, id = "gene")
  data.melt$value = log2(data.melt$value + 1)
  data.melt$ident = rep(cell.ident, each = length(gene))
  p = ggplot(data.melt, aes(factor(ident), value, fill = gene))
  p2 = p + geom_violin(scale = "width", adjust = adjust.use, trim = TRUE, aes(fill = factor(ident))) +
    ylab("Log2(cpm)")
  p2 = p2 + theme_bw()
  p3 = p2 + geom_jitter(height = 0, size = 0.1, colour = "purple", alpha = 0.1) +
    xlab("Cluster")
  p4 = p3 + theme(axis.title.x = element_text(face = "bold", colour = "#990000",
                                              size = 8), axis.text.x = element_text(vjust = 0.5, size = 8))
  p4 = p4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    guides(fill = FALSE)
  p5 = p4 + theme(axis.title.y = element_text(face = "bold", colour = "#990000",
                                              size = 8), axis.text.y = element_text(angle = 90, vjust = 0.5, size = 8))
  p5 = p5 + ggtitle(gene) + theme(plot.title = element_text(size = 12, face = "bold"))
  return(p5)
}


MultiPlotList <- function(plots, file, cols = 1, layout = NULL) {
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}


Multiple_Violin <- function(object, features.plot, ident.use, nCol = NULL) {
  if (is.null(nCol)) {
    nCol = 2
    if (length(features.plot) > 6)
      nCol = 3
    if (length(features.plot) > 9)
      nCol = 4
  }
  pList = lapply(features.plot, function(x) Violin_single(x, object[x, , drop = FALSE],
                                                          ident.use))
  MultiPlotList(pList, cols = nCol)
}

# To select a vector of genes present in the dataset

# Enter the gene names to be displayed in here

gene_selected_highlyExpressed_pluri <- c("SOX2", "POU5F1", "NANOG", "KLF4", "LEFTY2",
                                         "DNMT3B", "SDC2", "OTX2", "GDF3", "GFP42", "DPPA5", "DPPA2", "UTF1", "NODAL",
                                         "HHEX")

LASSO_1VS234 <- c("YBX1", "GNB2L1", "NGFRAP1", "HSPE1", "MCUR1", "FTH1", "PRELID1",
                  "HMGA1", "SERF2")

LASSO_4VS123 <- c("AURKAIP1", "NDUFS6", "SIVA1", "NMT1", "MAD2L2", "HSF2", "SEC11A",
                  "SLC25A1", "PTP4A2", "EZR", "PGP", "NRBP1", "FAM104B", "ACLY")

LASSO_3VS124 <- c("ZDHHC5", "TENM3", "NLGN4X", "COL4A2", "CASC4", "CABLES1", "SESN2",
                  "PLEC", "IRS2", "WRNIP1")

LASSO_2VS134 <- c("EIF1AY", "DYNC1LI1", "SRSF7", "GSDMD", "RNMTL1", "MAD2L1", "FKBP1A",
                  "FNTA", "NUDCD1", "SAPCD2")
gene_vec <- LASSO_4VS123

# check if the genes are in the dataset
nameV <- row.names(ori_dat_cpm)
gene_present <- unlist(sapply(gene_vec, function(x) {
  return(nameV[nameV == x])
}))

# plot the pannel here:

setwd("/Users/quan.nguyen/Documents/Powell_group_MacQuan/HiPSC/FullRun/PLOTs_violin_heatmap/")
pdf("Violin_LassoGenes_4vs123.pdf")
Multiple_Violin(ori_dat_cpm, gene_present, my.clusters)
dev.off()
