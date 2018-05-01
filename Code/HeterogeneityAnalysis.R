# This script performs heterogeneity analysis
# It uses expression matrix after sample-to-sample normalisation

library(edgeR)
library(reshape2)

# load clustering information
load("./my.clusters.Obj")
cluster1 <- which(my.clusters == 1)
cluster2 <- which(my.clusters == 2)
cluster3 <- which(my.clusters == 3)
cluster4 <- which(my.clusters == 4)
mtrx_exprs_C1 <- ori_dat_cpm[, cluster1]
mtrx_exprs_C2 <- ori_dat_cpm[, cluster2]
mtrx_exprs_C3 <- ori_dat_cpm[, cluster3]
mtrx_exprs_C4 <- ori_dat_cpm[, cluster4]

coef <- function(mtrx_exprs) {
  # remove zero rows and columns
  mtrx_exprs <- mtrx_exprs[apply(mtrx_exprs, 1, function(x) {
    any(x > 0)
  }), ]
  mtrx_exprs <- mtrx_exprs[, apply(mtrx_exprs, 2, function(x) {
    any(x > 0)
  })]

  # coefficient of variation calculated by standard deviation divided by the mean
  RowCoefVar <- function(x) {
    sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))/rowMeans(x)
  }

  variance_mtrx <- RowCoefVar(mtrx_exprs)

  variance_mtrx_MT <- as.data.frame(variance_mtrx)
  variance_mtrx_MT$GeneIDs <- row.names(mtrx_exprs)
  return(variance_mtrx_MT)
}

coeff_c1 <- coef(mtrx_exprs_C1)
colnames(coeff_c1) <- c("coeff_c1", "GeneIDs")
coeff_c2 <- coef(mtrx_exprs_C2)
colnames(coeff_c2) <- c("coeff_c2", "GeneIDs")
coeff_c3 <- coef(mtrx_exprs_C3)
colnames(coeff_c3) <- c("coeff_c3", "GeneIDs")
coeff_c4 <- coef(mtrx_exprs_C4)
colnames(coeff_c4) <- c("coeff_c4", "GeneIDs")
coeff_all <- coef(mtrx_exprs)
colnames(coeff_all) <- c("coeff_all", "GeneIDs")

merge12 <- merge(coeff_c1, coeff_c2, by.x = "GeneIDs", by.y = "GeneIDs")
merge34 <- merge(coeff_c3, coeff_c4, by.x = "GeneIDs", by.y = "GeneIDs")
merge1234 <- merge(merge12, merge34, by.x = "GeneIDs", by.y = "GeneIDs")
merge1234_all <- merge(merge1234, coeff_all, by.x = "GeneIDs", by.y = "GeneIDs")

colnames(merge1234_all) <- c("GeneIDs", "Cluster1", "Cluster2", "Cluster3", "Cluster4",
                             "AllCells")
merge1234_all <- merge1234_all[c(1, 6, 2:5)]

merge1234_all_melt <- melt(merge1234_all, vars.id = "GeneIDs")

size.x.use = 18
size.y.use = 16
size.title.use = 20
adjust.use = 1
size.use = 1

p = ggplot(merge1234_all_melt, aes(fill = factor(variable), y = value, x = factor(variable)))

p2 = p + geom_violin(scale = "width", adjust = adjust.use, trim = TRUE, aes(fill = factor(variable))) +
  ylab("Coefficient of variation")

p3 = p2 + guides(fill = FALSE) + geom_jitter(height = 0, size = 0.1, fill = "#990000",
                                             alpha = 0.1)
p4 = p3 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 = p4 + theme(axis.title.y = element_text(face = "bold", colour = "#990000", size = 16),
                axis.text.y = element_text(angle = 90, vjust = 0.5, size = 16)) + theme(axis.title.x = element_blank())
p5 = p5 + theme(axis.text.x = element_text(size = 12))
p5

# Cell cycle analysis-----------------------------------------------------------
library(scran)
# use normalised count data
load('dat_dcvl_mtrx_T_unLog_minus1_positive.Obj')
ddSeqdat <-t
# reformat names
names <-gsub(".*_", "", row.names(ddSeqdat))
# load reference
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
ddSeqExprsMatrix <-as.matrix(ddSeqdat)
CCassignments <- cyclone(ddSeqExprsMatrix, pairs=hs.pairs, names)
saveRDS(CCassignments, 'Cell_cycle_assignment.RDS')
# load my.clusters object
load('my.clusters.Obj')
# count cells
CC_clusters <- as.data.frame('phases'=CCassignments$phases, 'clusters'=my.clusters)
CC_clusters %>% count(clusters, phases)
CC_clusters %>% count(clusters, phases)
