#!/usr/bin/env python
# Core modules
import os

# Conda modules
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

# This sets the font used on the graph
font = {
    'family': 'Arial',
    'size': 10,
}

# Global stylings - set the font to above, white background, skinny lines
matplotlib.rc('font', **font)
matplotlib.rc("figure", facecolor="white")
matplotlib.rcParams['lines.linewidth'] = 0.5

def manage_subplot(subplot_obj, dataframe_dict):
    # Build the subplot
    build_series(dataframe_dict["Cluster 1"], "Cluster 1", "#e3ba22", subplot_obj)
    build_series(dataframe_dict["Cluster 2"], "Cluster 2", "#e6842a", subplot_obj)
    build_series(dataframe_dict["Cluster 3"], "Cluster 3", "#16687f", subplot_obj)
    build_series(dataframe_dict["Cluster 4"], "Cluster 4", "#8e6d8a", subplot_obj)

    # Set axis labels
    subplot_obj.set_xlabel('tSNE1')
    subplot_obj.set_ylabel('tSNE2')
    subplot_obj.set_zlabel('tSNE3')

    # Remove horrible grey background
    subplot_obj.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    subplot_obj.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    subplot_obj.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

def build_series(dataframe, label, colour, ax):
    x = dataframe['tSNE1']
    y = dataframe['tSNE2']
    z = dataframe['tSNE3']
    ax.scatter(x, y, z, c=colour, marker="o", lw=0.0, label=label, s=3)

# Import CSV data
data_path = "/Users/a.senabouth/Documents/IPSCPilot_scRNA/Publishing/data_file_figure1_3Dplot.txt"
info_df = pd.read_csv(data_path, sep="\t", index_col=0)
info_df['batchID'] = "Batch " + info_df['batchID'].astype(str)

# Convert into Dataframes
# Separate the datasets by clusterID
cluster1_data = info_df.loc[info_df['clusterID'] == 1]
cluster2_data = info_df.loc[info_df['clusterID'] == 2]
cluster3_data = info_df.loc[info_df['clusterID'] == 3]
cluster4_data = info_df.loc[info_df['clusterID'] == 4]

dataframe_dict = {"Cluster 1": cluster1_data, "Cluster 2": cluster2_data, "Cluster 3": cluster3_data, "Cluster 4": cluster4_data}
colour_dict = {"Cluster 1": "#e3ba22", "Cluster 2": "#e6842a", "Cluster 3": "#16687f", "Cluster 4": "#8e6d8a"}
fig = plt.figure(num=1, figsize=(11, 8.5), dpi=300)

# Normal View
ax1 = fig.add_subplot(221, projection='3d')
manage_subplot(ax1, dataframe_dict)

# XY-axis view
ax2 = fig.add_subplot(222, projection='3d')
manage_subplot(ax2, dataframe_dict)
ax2.view_init(0,90)
ax2.axes.get_yaxis().set_ticks([])
ax2.set_yticklabels([])

# XZ-axis view
ax3 = fig.add_subplot(223, projection='3d')
manage_subplot(ax3, dataframe_dict)
ax3.view_init(0,0)
ax3.axes.get_xaxis().set_ticks([])
ax3.set_xticklabels([])

# YZ-axis view
ax4 = fig.add_subplot(224, projection='3d')
manage_subplot(ax4, dataframe_dict)
ax4.view_init(90,0)
ax4.set_zticklabels([])

# Set legend
# We need to make "artist proxies" since figure legends are not supported for 3D subplots :()
cluster1_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#e3ba22")
cluster2_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#e6842a")
cluster3_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#16687f")
cluster4_proxy = plt.Line2D([0], [0], linestyle="none", marker="o", markersize=10, markerfacecolor="#8e6d8a")
fig.legend((cluster1_proxy, cluster2_proxy, cluster3_proxy, cluster4_proxy), ("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), loc="lower center", numpoints=1, ncol=4)
fig.tight_layout()
# fig.tight_layout(pad=0.2, w_pad=0.4, h_pad=0.8)
#Create a pdf file to save to.
pp = PdfPages("Figure1a.pdf")
plt.savefig(pp, format="pdf")
# plt.show()

pp.close()
