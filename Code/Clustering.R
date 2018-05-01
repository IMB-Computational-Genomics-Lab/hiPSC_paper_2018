# This cript runs the main clustering pipeline.  This script was later developed and
# implemented in the ascend R package, available at https://github.com/IMB-Computational-Genomics-Lab/ascend
# This script outputs the most stable cluster the first 25% search space

################################################
# Running 10 clustering windows, covering 25% of the total search space to find
# stable clustering result
################################################

library(dynamicTreeCut)
# Load the input data, which is the PCA transformed data (using the first 10 PCs)
load("Matrix_PCAreduced_forHC_clust.Obj")
sce_dcvl_PCA_transformed <- PCA_customPL
dist_t_exprs <- dist(sce_dcvl_PCA_transformed)
saveRDS(dist_t_exprs, "shared_distance_matrix.RDS")

# The original tree is the shared tree to find the optimal clusters

original.tree <- hclust(dist_t_exprs, method = "ward.D2")
# The original clusters to be used as the reference
original.clusters <- unname(cutreeDynamic(original.tree, distM = as.matrix(dist_t_exprs),
                                          verbose = 0))

saveRDS(original.clusters, "my.clusters_defaul_no_minSplitHeight.RDS")
table(original.clusters)
original.tree$labels <- rep("", length(original.clusters))
plot(original.tree)

# A loop for running clustering using a range of different tree cutting parameters
# Results are saved into the working folder

clustering_param <- function(Height = 0) {
  my.clusters <- unname(cutreeDynamic(original.tree, distM = as.matrix(dist_t_exprs),
                                      minSplitHeight = Height, verbose = 0))
  saveRDS(my.clusters, paste0("my.clusters_", Height, ".RDS"))
}

# run 10 clustering rounds
for (i in seq(0.025:0.25, by = 0.025)) {
  clustering_param(i)
}

################################################
# Plotting the dendrogram with clustering results from running 10 different tree
# height parameter values
################################################

library(dynamicTreeCut)
dist_t_exprs <- readRDS("shared_distance_matrix.RDS")

# The original tree to be used as the reference
original.tree <- hclust(dist_t_exprs, method = "ward.D2")
# The original clusters to be used as the reference
original.clusters <- unname(cutreeDynamic(original.tree, distM = as.matrix(dist_t_exprs),
                                          verbose = 0))
saveRDS(original.clusters, "my.clusters_defaul_no_minSplitHeight.RDS")

# load the saved clustering results from the above 10 clustering runs
list_clusters <- list()
for (i in seq(0.025:0.25, by = 0.025)) {
  list_clusters <- c(list_clusters, paste0("my.clusters_", i, ".RDS"))
}

# read the first 10 clusters
read_clusters <- list()
for (i in 1:10) {
  read_clusters[[i]] <- readRDS(list_clusters[[i]])
}

col_all <- matrix(unlist(read_clusters), ncol = 10)
col_all <- as.data.frame(col_all)

read_clusters_ref <- readRDS("my.clusters_defaul_no_minSplitHeight.RDS")
col_all$ref <- read_clusters_ref

colnames(col_all) <- c(seq(1:10), "ref")

original.tree$labels <- rep("", length(original.tree$labels))

# plot tree with clustering result bars underneath


#-----------------------------------------------------------------------------
# Function to plot dendrogram and color
#-----------------------------------------------------------------------------
plot_CORE2 <- function(original.tree, list_clusters = NULL, color_branch = NULL) {
  library(RColorBrewer)
  n <- length(unique(unlist(list_clusters[[1]])))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual", ]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (is.null(color_branch)) {
    color_branch <- col_vector
  }

  col_all <- matrix(unlist(list_clusters), ncol = length(list_clusters))
  col_all <- as.data.frame(col_all)
  # remove branch labels
  original.tree$labels <- rep("", length(original.tree$labels))

  plotDendroAndColors <- function(dendro, colors, groupLabels = NULL, rowText = NULL,
                                  rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, textPositions = NULL,
                                  setLayout = TRUE, autoColorHeight = TRUE, colorHeight = 0.2, rowWidths = NULL,
                                  dendroLabels = NULL, addGuide = FALSE, guideAll = FALSE, guideCount = 50,
                                  guideHang = 0.2, addTextGuide = FALSE, cex.colorLabels = 0.8, cex.dendroLabels = 0.9,
                                  cex.rowText = 0.8, marAll = c(1, 5, 3, 1), saveMar = TRUE, abHeight = NULL,
                                  abCol = "red", ...) {
    oldMar = par("mar")
    if (!is.null(dim(colors))) {
      nRows = dim(colors)[2]
    } else nRows = 1
    if (!is.null(rowText))
      nRows = nRows + if (is.null(textPositions))
        nRows else length(textPositions)
    if (autoColorHeight)
      colorHeight = 0.2 + 0.3 * (1 - exp(-(nRows - 1)/6))
    if (setLayout)
      layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight, colorHeight))
    par(mar = c(0, marAll[2], marAll[3], marAll[4]))

    # quan add this for a case--------------------
    dendro1 <- as.dendrogram(dendro)
    dendro2 <- dendro1 %>% set("branches_k_color", k = 4, value = c("#916F8E",
                                                                    "#E48637", "#216B7F", "#E3B939"))
    plot(dendro2, labels = dendroLabels, cex = cex.dendroLabels, ...)
    # done the addition---------------------------

    if (addGuide)
      addGuideLines(dendro, count = if (guideAll)
        length(dendro$height) + 1 else guideCount, hang = guideHang)
    if (!is.null(abHeight))
      abline(h = abHeight, col = abCol)
    par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
    plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels,
                       rowText = rowText, rowTextAlignment = rowTextAlignment, rowTextIgnore = rowTextIgnore,
                       textPositions = textPositions, cex.rowText = cex.rowText, rowWidths = rowWidths,
                       addTextGuide = addTextGuide)
    if (saveMar)
      par(mar = oldMar)
  }
  #-----------------------------------------------------------------------------
  # plotColorUnderTree
  #-----------------------------------------------------------------------------
  plotColorUnderTree <- function(dendro, colors, rowLabels = NULL, rowWidths = NULL,
                                 rowText = NULL, rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL,
                                 textPositions = NULL, addTextGuide = TRUE, cex.rowLabels = 1, cex.rowText = 0.8,
                                 ...) {

    plotOrderedColors(dendro$order, colors = colors, rowLabels = rowLabels, rowWidths = rowWidths,
                      rowText = rowText, rowTextAlignment = rowTextAlignment, rowTextIgnore = rowTextIgnore,
                      textPositions = textPositions, addTextGuide = addTextGuide, cex.rowLabels = cex.rowLabels,
                      cex.rowText = cex.rowText, startAt = 0, ...)
  }
  #-----------------------------------------------------------------------------
  # plotOrderedColors
  #-----------------------------------------------------------------------------
  plotOrderedColors <- function(order, colors, rowLabels = NULL, rowWidths = NULL,
                                rowText = NULL, rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL,
                                textPositions = NULL, addTextGuide = TRUE, cex.rowLabels = 1, cex.rowText = 0.8,
                                startAt = 0, ...) {
    colors = as.matrix(colors)
    dimC = dim(colors)
    if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) == dimC[2]))
      rowLabels = colnames(colors)
    sAF = options("stringsAsFactors")
    options(stringsAsFactors = FALSE)
    on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)
    nColorRows = dimC[2]
    if (length(order) != dimC[1])
      stop("ERROR: length of colors vector not compatible with number of objects in 'order'.")
    C = colors[order, , drop = FALSE]
    step = 1/(dimC[1] - 1 + 2 * startAt)
    barplot(height = 1, col = "white", border = FALSE, space = 0, axes = FALSE)
    charWidth = strwidth("W")/2
    if (!is.null(rowText)) {
      if (is.null(textPositions))
        textPositions = c(1:nColorRows)
      if (is.logical(textPositions))
        textPositions = c(1:nColorRows)[textPositions]
      nTextRows = length(textPositions)
    } else nTextRows = 0
    nRows = nColorRows + nTextRows
    ystep = 1/nRows
    if (is.null(rowWidths)) {
      rowWidths = rep(ystep, nColorRows + nTextRows)
    } else {
      if (length(rowWidths) != nRows)
        stop("plotOrderedColors: Length of 'rowWidths' must equal the total number of rows.")
      rowWidths = rowWidths/sum(rowWidths)
    }
    hasText = rep(0, nColorRows)
    hasText[textPositions] = 1
    csPosition = cumsum(c(0, hasText[-nColorRows]))
    colorRows = c(1:nColorRows) + csPosition
    rowType = rep(2, nRows)
    rowType[colorRows] = 1
    physicalTextRow = c(1:nRows)[rowType == 2]
    yBottom = c(0, cumsum(rowWidths[nRows:1]))
    yTop = cumsum(rowWidths[nRows:1])
    if (!is.null(rowText)) {
      rowTextAlignment = match.arg(rowTextAlignment)
      rowText = as.matrix(rowText)
      textPos = list()
      textPosY = list()
      textLevs = list()
      for (tr in 1:nTextRows) {
        charHeight = max(strheight(rowText[, tr], cex = cex.rowText))
        width1 = rowWidths[physicalTextRow[tr]]
        nCharFit = floor(width1/charHeight/1.7/par("lheight"))
        if (nCharFit < 1)
          stop("Rows are too narrow to fit text. Consider decreasing cex.rowText.")
        set = textPositions[tr]
        textLevs[[tr]] = sort(unique(rowText[, tr]))
        textLevs[[tr]] = textLevs[[tr]][!textLevs[[tr]] %in% rowTextIgnore]
        nLevs = length(textLevs[[tr]])
        textPos[[tr]] = rep(0, nLevs)
        orderedText = rowText[order, tr]
        for (cl in 1:nLevs) {
          ind = orderedText == textLevs[[tr]][cl]
          sind = ind[-1]
          ind1 = ind[-length(ind)]
          starts = c(if (ind[1]) 1 else NULL, which(!ind1 & sind) + 1)
          ends = which(c(ind1 & !sind, ind[length(ind)]))
          if (length(starts) == 0)
            starts = 1
          if (length(ends) == 0)
            ends = length(ind)
          if (ends[1] < starts[1])
            starts = c(1, starts)
          if (ends[length(ends)] < starts[length(starts)])
            ends = c(ends, length(ind))
          lengths = ends - starts
          long = which.max(lengths)
          textPos[[tr]][cl] = switch(rowTextAlignment, left = starts[long],
                                     center = (starts[long] + ends[long])/2 + 0.5, right = ends[long] +
                                       1)
        }
        if (rowTextAlignment == "left") {
          yPos = seq(from = 1, to = nCharFit, by = 1)/(nCharFit + 1)
        } else {
          yPos = seq(from = nCharFit, to = 1, by = -1)/(nCharFit + 1)
        }
        textPosY[[tr]] = rep(yPos, ceiling(nLevs/nCharFit) + 5)[1:nLevs][rank(textPos[[tr]])]
      }
    }
    jIndex = nRows
    if (is.null(rowLabels))
      rowLabels = c(1:nColorRows)
    C[is.na(C)] = "grey"
    for (j in 1:nColorRows) {
      jj = jIndex
      ind = (1:dimC[1])
      xl = (ind - 1.5 + startAt) * step
      xr = (ind - 0.5 + startAt) * step
      yb = rep(yBottom[jj], dimC[1])
      yt = rep(yTop[jj], dimC[1])
      if (is.null(dim(C))) {
        rect(xl, yb, xr, yt, col = as.character(C), border = as.character(C))
      } else {
        rect(xl, yb, xr, yt, col = as.character(C[, j]), border = as.character(C[,
                                                                                 j]))
      }
      text(rowLabels[j], pos = 2, x = -charWidth/2 + xl[1], y = (yBottom[jj] +
                              yTop[jj])/2, cex = cex.rowLabels, xpd = TRUE)
      textRow = match(j, textPositions)
      if (is.finite(textRow)) {
        jIndex = jIndex - 1
        xt = (textPos[[textRow]] - 1.5) * step
        xt[xt < par("usr")[1]] = par("usr")[1]
        xt[xt > par("usr")[2]] = par("usr")[2]
        yt = yBottom[jIndex] + (yTop[jIndex] - yBottom[jIndex]) * (textPosY[[textRow]] +
                                                                     1/(2 * nCharFit + 2))
        nt = length(textLevs[[textRow]])
        if (addTextGuide)
          for (l in 1:nt) lines(c(xt[l], xt[l]), c(yt[l], yTop[jIndex]), col = "darkgrey",
                                lty = 3)
        textAdj = c(0, 0.5, 1)[match(rowTextAlignment, c("left", "center",
                                                         "right"))]
        text(textLevs[[textRow]], x = xt, y = yt, adj = c(textAdj, 1), xpd = TRUE,
             cex = cex.rowText)
      }
      jIndex = jIndex - 1
    }
    for (j in 0:(nColorRows + nTextRows)) lines(x = c(0, 1), y = c(yBottom[j +
         1], yBottom[j + 1]))
  }

  #-----------------------------------------------------------------------------
  # Start selecting colors
  #-----------------------------------------------------------------------------

  color_range <- color_branch
  number_colors <- length(unique(unlist(col_all)))
  col_all2 <- col_all
  for (i in 1:number_colors) {
    col_all2[col_all2 == i] <- as.character(color_range[i])
  }
  colnames(col_all2) <- gsub("V", "", colnames(col_all2))

  plotDendroAndColors(original.tree, col_all2)

}

# run the plotCORE2 function to generate the supplementary figure
pdf("CORE_clustering_10windows_SupplementaryFigure.pdf")
plot_CORE2(original.tree = original.tree, list_clusters = read_clusters[1:10],
           color_branch = c("#E3B939",  "#E48637", "#216B7F", "#916F8E"))
dev.off()
