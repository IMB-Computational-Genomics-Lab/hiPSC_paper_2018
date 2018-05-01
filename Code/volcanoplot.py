#!/usr/bin/env python
import os

# Getting around a matplotlib issue - using TkAgg backend instead of MacOSX
import matplotlib
matplotlib.use('TkAgg')

# Anaconda packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from adjustText import adjust_text
from matplotlib.backends.backend_pdf import PdfPages

####Known marker genes#######
PL_list = set(["CXCL5", "IDO1", "LCK", "TRIM22", "DNMT3B", "HESX1", "SOX2", "POU5F1", "NANOG", 'TBX3', 'KLF4', 'KLF5'])
EC_list = set(["POU4F1", "TRPM8", "OLFM3", "DMBX1", "CDH9", "NOS2", "MYO3B", "PAPLN", "DRD4", "WNT1", "LMX1A","ZBTB16", "NR2F1/2", "PAX3", "PAX6","MAP2","COL2A1","SOX1","NR2F2", "SDC2","EN1"])
ME_list = set(["ODAM", "ESM1", "HOPX", "PLVAP", "FOXF1", "TM4SF1", "FCN3", "COLEC10","HAND2", "HAND1", "CDX2","RGS4", "BMP10", "TBX3", "SST", "ABCA4","NKX2-5", "HEY1", "PDGFRA", "SNAI2", "KLF5", "GATA4", "CDH5", "IL6ST", "ALOX5", "MIXL1", 'POU5F1', 'T'])
EN_list = set(["FOXP2", "HMP19", "CDH20", "ELAVL3", "PHOX2B", "CPLX2", "POU3F", "CABP7", "LEFTY1", "EOMES", "NODAL", "LEFTY2", "FOXA2", "HNF4A", "RXRG", "HNF1B", "SOX17", "HHEX", "GATA6", "PRDM1", "FOXA1", "CLDN1"])
MS_list = set(["T", "FGF4", "GDF3", "NR5A2", "PTHLH", "NPPB"])

ALL_KNOWNMARKERS = list(set().union(PL_list, EC_list, ME_list, EN_list, MS_list))

# Custom labelling algorithm
def labelpoints_custom(dataframe_obj, subplot_obj, label_df):
    # Iterate through label_df to determine if they get a label or not
    point_list = {}

    for idx, booleon_tup in label_df.iterrows():
        # Examine their labels
        x_bool, y_bool = booleon_tup.tolist()
        # Grab the x,y coordinates of the point
        gene, x, y = (dataframe_obj.loc[idx][['id', 'log2FoldChange', 'padj']]).tolist()
        # If the gene is in known markers and are outliers
        if gene in ALL_KNOWNMARKERS:
            if y > 5:
                point_list[gene] = (x, y)
        else:
            group = sortgroup(x_bool, y_bool)
            if group in (3, 4):
                point_list[gene] = (x, y)


def labelpoints_module(dataframe_obj, subplot_obj, label_df):
    texts = []
    for idx, booleon_tup in label_df.iterrows():
        x_bool, y_bool = booleon_tup.tolist()
        gene, x, y = (dataframe_obj.loc[idx][['id', 'log2FoldChange', 'padj']]).tolist()
        if gene in ALL_KNOWNMARKERS:
            if y > 5:
                text_obj = subplot_obj.text(x, y, gene, size=6)
                text_obj.set_path_effects([path_effects.Stroke(linewidth=1.8, foreground="white"), path_effects.Normal()])
                texts.append(text_obj)
        else:
            group = sortgroup(x_bool, y_bool)
            if group in (3, 4):
                text_obj = subplot_obj.text(x, y, gene, size=6)
                text_obj.set_path_effects([path_effects.Stroke(linewidth=1.8, foreground="white"), path_effects.Normal()])
                texts.append(text_obj)

    adjust_text(texts, precision=0, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

# Finding limits to set thresholds for coloured and labelled groups
def findlimits(x_series, y_series, cutoff):
    x_lim1, x_lim2 = x_series.quantile([cutoff, 1-cutoff])
    y_lim1, y_lim2 = y_series.quantile([cutoff, 1-cutoff])
    if y_lim1 == (-0.0):
        y_lim1 = 0.0

    return (x_lim1, x_lim2), (y_lim1, y_lim2)

# Point annotation using adjust_text module
def annotatePoints(dataframe_obj, subplot_obj):
    x_limits, y_limits = findlimits(dataframe_obj['log2FoldChange'], dataframe_obj['padj'], 0.0005)
    x_sorted = dataframe_obj['log2FoldChange'].map(lambda x: testboundaries(x_limits, x))
    y_sorted = dataframe_obj['padj'].map(lambda y: testboundaries(y_limits, y))

    label_df = pd.concat([x_sorted, y_sorted], axis=1)
    # Test and compile list of points that should be labelled
    labelpoints_module(dataframe_obj, subplot_obj, label_df)

def buildVolcano(dataframe_obj, subplot_obj, adjust=False):
    # Extract sub-datasets
    # Group 1 - normal samples
    group1_df = dataframe_obj[dataframe_obj["group"] == 1]
    # Group 2 - Known markers
    group2_df = dataframe_obj[dataframe_obj["group"] == 2]
    # Group 3 - Significant along X axis
    group3_df = dataframe_obj[dataframe_obj["group"] == 3]
    # Group 4 - Significant along Y axis
    group4_df = dataframe_obj[dataframe_obj["group"] == 4]

    # Plot
    # subplot_obj.scatter(dataframe_obj['log2FoldChange'], dataframe_obj['padj'], label='Normal', color="black", s=10)
    subplot_obj.scatter(group1_df['log2FoldChange'], group1_df['padj'], label='Normal', color="#b2b2b2", s=8)
    subplot_obj.scatter(group2_df['log2FoldChange'], group2_df['padj'], label='Known markers', color="#ab775e", s=8)
    subplot_obj.scatter(group3_df['log2FoldChange'], group3_df['padj'], label='Top Fold-Change', color="#673f5c", s=8)
    subplot_obj.scatter(group4_df['log2FoldChange'], group4_df['padj'], label='Top P-Values', color="#394e64", s=8)

def sortgroup(x_bool, y_bool):
    bool_obj = (x_bool, y_bool)

    if True not in bool_obj:
        return 1
    else:
        if y_bool:
            if x_bool:
                return 3
            else:
                return 4
        else:
            return 3

def testboundaries(limit_tup, value):
    if not (limit_tup[0] < value < limit_tup[1]):
        if limit_tup[0] == value:
            return False
        else:
            return True
    else:
        return False

def buildSubplot(figure_obj, subplot_obj, dataframe_obj, title):
    subplot_obj.set_xlabel("log2FoldChange", fontsize=8, fontweight='bold')
    subplot_obj.set_ylabel("-log10(P-value)", fontsize=8, fontweight='bold')
    subplot_obj.spines['right'].set_visible(False)
    subplot_obj.spines['top'].set_visible(False)
    subplot_obj.xaxis.set_ticks_position('bottom')
    subplot_obj.yaxis.set_ticks_position('left')
    subplot_obj.tick_params(axis="both", which="major", labelsize=8)
    subplot_obj.tick_params(axis="both", which="minor", labelsize=8)
    subplot_obj.set_title(title, fontsize=9, fontweight='bold')
    buildVolcano(dataframe_obj, subplot_obj)
    annotatePoints(dataframe_obj, subplot_obj)

def parseDatasets(dataframe_obj):
    # Convert Pvals into log10
    dataframe_obj['padj'] = -(np.log(dataframe_obj['padj']))
    dataframe_obj = dataframe_obj.replace(to_replace=-0.000000, value=0.000000)
    dataframe_obj = dataframe_obj.dropna()

    # Define limits on what is significant
    x_limits, y_limits = findlimits(dataframe_obj['log2FoldChange'], dataframe_obj['padj'], 0.005)
    # Create a new column
    series_column = pd.Series()
    # Group datasets
    for idx, row in dataframe_obj.iterrows():
       if row['id'] in ALL_KNOWNMARKERS:
           series_column = series_column.set_value(idx, 2)
       else:
           x_booleon = testboundaries(x_limits, row['log2FoldChange'])
           y_booleon = testboundaries(y_limits, row['padj'])
           group_val = sortgroup(x_booleon, y_booleon)
           series_column = series_column.set_value(idx, group_val)
    dataframe_obj['group'] = series_column
    return dataframe_obj

# Pandas Import Function
def importcsv(path):
    dataset = pd.read_csv(path, delimiter='\t')
    processed_dataset = parseDatasets(dataset)
    return processed_dataset

# Filepaths for input files
working_dir = "/Users/a.senabouth/Documents/IPSCPilot_scRNA/ProcessedData/Diagram Generation/Figure1c/"
cluster1_dir = working_dir + "DEgenes_for_VolcanoPlot_plotly_C1vsC234.txt"
cluster2_dir = working_dir + "DEgenes_for_VolcanoPlot_plotly_C2vsC124.txt"
cluster3_dir = working_dir + "DEgenes_for_VolcanoPlot_plotly_C3vsC124.txt"
cluster4_dir = working_dir + "DEgenes_for_VolcanoPlot_plotly_C4vsC123.txt"

# Plot figure
fig = plt.figure(num=1, figsize=(6, 6), dpi=600)

# Cluster 1 vs All
cluster1_df = importcsv(cluster1_dir)
ax1 = fig.add_subplot(221)
buildSubplot(fig, ax1, cluster1_df, "DE genes Cluster 1 vs Clusters 2,3 and 4")

# Cluster 2 vs All
cluster2_df = importcsv(cluster2_dir)
ax2 = fig.add_subplot(222)
buildSubplot(fig, ax2, cluster2_df, "DE genes Cluster 2 vs Clusters 1,3 and 4")

# Cluster 3 vs All
cluster3_df = importcsv(cluster3_dir)
ax3 = fig.add_subplot(223)
buildSubplot(fig, ax3, cluster3_df, "DE genes Cluster 3 vs Clusters 1,2 and 4")

# Cluster 4 vs All
cluster4_df = importcsv(cluster4_dir)
ax4 = fig.add_subplot(224)
buildSubplot(fig, ax4, cluster4_df, "DE genes Cluster 4 vs Clusters 1,2 and 3")

# Proxies for legend
normal_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#000000")
known_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#ff9538")
fold_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#930052")
pval_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#013b83")
#fig.legend((normal_proxy, known_proxy, fold_proxy, pval_proxy), ("Normal", "Known Markers", "Top Fold-Change", "Top P-Values"), loc="center", numpoints=1, ncol=4)

# Formatting for the entire figure
#fig.tight_layout()
fig.tight_layout()
# plt.show()
pp = PdfPages("Figure1c.pdf")
plt.savefig(pp, format="pdf")
pp.close()
