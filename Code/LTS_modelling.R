# This script is to run the LTS prediction pipeline.This cript is to estimate the best
# multigenic model for one cluster, compared to all the remaining clusters.
# The algorithm was later developed and implemented in the scGPS R package - available at https://github.com/IMB-Computational-Genomics-Lab/scGPS).
# For each prediction, it run 100 bootstraps.
library(glmnet)
library("foreach")
library("doParallel")
n_cores = 6

# register the cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

argv <- commandArgs(TRUE)

c_selectID = argv[1]
DE_result_file_ID = argv[2]

# load expression matrix (ori_dat)
load("dat_dcvl_mtrx_T_unLog_minus1_positive_BOC_for_edgeR.Obj")
ori_dat <- t

# load my.clusters object
load("my.clusters.Obj")

cat(unique(my.clusters))

names <- rownames(ori_dat)

names <- gsub("_.*", "", names)

rownames(ori_dat) <- names
# get cluster IDs
n_clusters <- length(unique(my.clusters))
# extract clusters

cluster_select <- which(my.clusters == as.numeric(c_selectID))

cluster_compare <- which(my.clusters != as.numeric(c_selectID))

cellNames_cluster <- cbind(colnames(ori_dat), my.clusters)

# to extract DE genes
DE_result <- read.table(paste0("result_table_DESeq_", c_selectID, "_vs_", DE_result_file_ID,
                               "_significant10x15K.txt"), header = T)

DEgenes <- DE_result$id
DEgenes <- gsub("_.*", "", DEgenes)
DE_idx <- which(rownames(ori_dat) %in% DEgenes)

Lit_New_Lasso <- function(cluster_select, c_selectID = "1", cluster_compare, c_compareID = "2",
                          M_or_DE_idx) {
  cluster_select_indx_S1 <- sample(cluster_select, round(length(cluster_select)/2),
                                   replace = F)
  # taking 50% out for training
  cluster_compare_indx_S1 <- sample(cluster_compare, round(length(cluster_compare)/2),
                                    replace = F)
  # prepare predictor matrix containing both clutering classes
  genes_S1 <- ori_dat[M_or_DE_idx, c(cluster_select_indx_S1, cluster_compare_indx_S1)]
  ### Selecting gene markers####
  predictor_S1 <- genes_S1
  # generate categorical class response set all values to cluster select (character
  # type)
  y_cat = rep(c_selectID, length(predictor_S1[1, ]))
  # replace values for cluster compare
  ori_compare <- ori_dat[, cluster_compare]
  # no need to change predictor_markerS1 here
  Sub_clustercompare_Indx_S1 <- which(colnames(predictor_S1) %in% colnames(ori_compare))
  # change value of the cluster id
  y_cat[Sub_clustercompare_Indx_S1] <- rep(c_compareID, length(Sub_clustercompare_Indx_S1))
  # tranpose prediction matrix
  predictor_S1 <- t(predictor_S1)
  # fit the lasso model
  fit <- glmnet(predictor_S1, y_cat, family = "binomial")
  # Prepare validation test; keep all cells except for those used in the training
  # set
  cluster_select_indx_S2 <- sample(cluster_select[-cluster_select_indx_S1])
  cluster_compare_indx_S2 <- sample(cluster_compare[-cluster_compare_indx_S1])
  genes_S2 <- ori_dat[M_or_DE_idx, c(cluster_select_indx_S2, cluster_compare_indx_S2)]
  predictor_S2 <- t(genes_S2)
  # fitting with cross validation
  cvfit = cv.glmnet(t(genes_S1), y_cat, family = "binomial", type.measure = "class")
  predict_clusters <- predict(cvfit, newx = predictor_S2, type = "class", s = cvfit$lambda.min)
  # to extract coefficient Beta for a gene for an optimized lambda value
  cvfit_out <- as.matrix(coef(cvfit, s = cvfit$lambda.min))
  cvfit_out <- as.data.frame(cvfit_out)
  # find number of genes with coefficient different to 0
  cvfit_out$name <- row.names(cvfit_out)
  sub_cvfit_out <- cvfit_out[cvfit_out$`1` != 0, ]

  # Extracting deviance explained
  t_DE <- as.matrix(print(cvfit$glmnet.fit))
  dat_DE <- as.data.frame(t_DE)
  colnames(dat_DE) <- c("Dfd", "Deviance", "lambda")
  # to get the coordinate for lambda that produces minimum error
  dat_DE_Lambda_idx <- which(round(dat_DE$lambda, digit = 3) == round(cvfit$lambda.min,
                                                                      digits = 3))
  dat_DE <- dat_DE[1:dat_DE_Lambda_idx[1], ]
  require(dplyr)
  dat_DE_fm_DE <- dat_DE %>% group_by(Dfd) %>% summarise(Deviance = max(Deviance))
  dat_DE_fm_DE <- as.data.frame(dat_DE_fm_DE)
  dat_DE_fm_DE$DEgenes <- paste0("DEgenes_C", c_selectID, "_day_0")
  remaining <- c("remaining", 1, "DEgenes")
  dat_DE_fm_DE <- rbind(dat_DE_fm_DE, remaining)
  # to return the output as 5 lists
  return(list(predict_clusters, predictor_S2, sub_cvfit_out, dat_DE_fm_DE, cvfit))
}

# a function to combine results from parallel runs
comb <- function(x, ...) {
  lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# running parrallel prediction with 100 bootstraps for DE genes
All_List <- foreach(i = 1:100, .combine = comb, .multicombine = TRUE, .init = list(list(),
          list(), list(), list()), .packages = "glmnet") %dopar% {
            # set the character values (one of the two classes) for the response variable
            c_compareID = 1:length(unique(my.clusters))
            c_compareID <- paste0(c_compareID[-which(c_compareID == c_selectID)], collapse = "")
            predict_marker <- Lit_New_Lasso(cluster_select, c_selectID = c_selectID, cluster_compare,
            c_compareID, M_or_DE_idx = DE_idx)
            predictor_S2_name <- predict_marker[[2]]
            predict_label <- predict_marker[[1]]
            # from here compare to original clusering classes to check for accuracy
            predict_index <- which(cellNames_cluster[, 1] %in% row.names(predictor_S2_name))
            original_cluster <- cellNames_cluster[predict_index, ]
            original_cluster <- original_cluster[order(original_cluster[, 2], decreasing = T),
            ]
            original_cluster <- as.data.frame(original_cluster)
            predict_clusters <- as.data.frame(predict_label)
            predict_clusters$cellnames <- row.names(predict_clusters)
            compare <- merge(original_cluster, predict_clusters, by.x = "V1", by.y = "cellnames")
            # change cluster IDs here
            cluster_select_predict <- subset(compare, (as.numeric(compare$my.clusters) ==
                                             c_selectID & compare$`1` == c_selectID) | (compare$my.clusters != c_selectID &
                                             compare$`1` != c_selectID))
            accurate <- dim(cluster_select_predict)[1]
            inaccurate <- dim(compare)[1] - dim(cluster_select_predict)[1]
            list_acc_inacc <- list(accurate, inaccurate)
            # Return the list for the combine function during the parallelisation process
            return(list(list_acc_inacc = list_acc_inacc, list_SigGenes = predict_marker[[3]],
                  list_Deviance = predict_marker[[4]], list_cvFit = predict_marker[[5]]))
                  saveRDS(i, paste0("Temp_reporting_Done_BootStrap_", i))
          }


saveRDS(All_List, file = paste0("All_List_LASSO_for_Cluster", c_selectID, "vs_", DE_result_file_ID,
                                ".RDS"))
cat("Done Boostrap DEgenes")

stopCluster(cl)

# Summarise bootstrap results and perform prediction############################

library(glmnet)
library(dplyr)
list.files()

# load data from the bootstrap output
dat <- readRDS("All_List_LASSO_for_Cluster3vs_cluster234.RDS")
dat <- readRDS("All_List_LASSO_for_Cluster2vs_cluster134.RDS")
dat <- readRDS("All_List_LASSO_for_Cluster1vs_cluster234.RDS")
dat <- readRDS("All_List_LASSO_for_Cluster4vs_cluster123.RDS")

inac_accr <- dat[[1]]
genesSig <- dat[[2]]
deviDat <- dat[[3]]
cvfit100 <- dat[[4]]

percent_inacc <- c()
for (i in 1:length(inac_accr)) {
  percent_inacc[i] <- as.numeric(inac_accr[[i]][2])/(as.numeric(inac_accr[[i]][2]) +
                                                       as.numeric(inac_accr[[i]][1])) * 100
}
percent_inacc <- as.data.frame(percent_inacc)
percent_inacc$order <- row.names(percent_inacc)

# use deviDat and cvfit100 to check order of the paralleled outputs
percent_inacc
deviance <- vector()
for (i in 1:length(deviDat)) {
  deviance <- c(deviance, deviDat[[i]]$Deviance[dim(deviDat[[i]])[1] - 1])
}
acc_inacc_dev <- percent_inacc
acc_inacc_dev$Deviance <- deviance
acc_inacc_dev

library(dplyr)
acc_inacc_dev_od <- acc_inacc_dev %>% arrange(percent_inacc, Deviance)

# can plot pecent inaccurate vs deviance explained to see the linear relationship
best_bootstrap <- as.numeric(acc_inacc_dev_od$order[1])
genesSig[[best_bootstrap]]
head(acc_inacc_dev_od)

# perform prediction

# select a timepoint
dayID1 = "day0"

cvfit_best <- cvfit100[[best_bootstrap]]
class(cvfit_best) <- c("cv.glmnet")
cvfit_best
genesSig[[best_bootstrap]]
# saveRDS(cvfit_best, paste0('cvfit_best_for_cluster_',c_selectID_1,
# '_onDay_',dayID1))

# load clusters
load("../my.clusters.Obj")
table(my.clusters)
# load expression matrix
ori_dat_2 <- readRDS("../Expression_CPM_unLog_minus1_positive.Obj")

c_selectID_2 = 1

cluster_select <- which(my.clusters == as.numeric(c_selectID_2))
length(cluster_select)
names <- gsub("_.*", "", rownames(ori_dat_2))
rownames(ori_dat_2) <- names
head(names)
gene_cvfit <- cvfit_best$glmnet.fit$beta@Dimnames[[1]]
gene_cvfit
cvfitGenes_idx <- which(rownames(ori_dat_2) %in% gene_cvfit)
length(gene_cvfit)
length(cvfitGenes_idx)
to_add <- length(cvfitGenes_idx) - length(gene_cvfit)
to_add_idx <- sample(cvfitGenes_idx, to_add, replace = F)
predictor_S2_temp <- ori_dat_2[c(to_add_idx, cvfitGenes_idx), cluster_select]
predict_clusters <- predict(cvfit_best, newx = t(predictor_S2_temp), type = "class",
                            s = cvfit_best$lambda.min)
predict_clusters <- as.data.frame(predict_clusters)
predict_clusters$names <- gsub("-.$", "", row.names(predict_clusters))
predict_clusters$names <- gsub("^._", "", predict_clusters$names)
colnames(predict_clusters) <- c("predicted_cluster", "names")
head(predict_clusters)

# get the results as the percent of cells belonging to each of the two classes
table(predict_clusters$predicted_cluster)
