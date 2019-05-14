import pandas as pd
import numpy as np
from scipy.spatial import distance
import networkx as nx
import string
import plotly.plotly as py
from plotly.graph_objs import *
import plotly
import matplotlib.pyplot as plt
import seaborn as sns

plotly.tools.set_credentials_file(username='slt666666', api_key='KEROY4DT9JraomPybSjQ')

# # extract only "gene" list from gff
# gff_tomato = pd.read_table("original_data/ITAG3.2_gene_models.gff", header=None, skiprows=14)
# gff_tomato = gff_tomato[gff_tomato.iloc[:, 2] == "gene"]
# # print(gff_tomato.head())
#
# # extract NLR gene_ids from data_AKI
# NLR_list = pd.read_csv("original_data/tomato_NLRs_from_AKI.csv")
# NLR_list = np.sort(NLR_list.gene_short_na)
# # print(NLR_list[0:10])
#
# # extract gene info
# NLR_gene_info = pd.DataFrame()
# for gene_id in NLR_list:
#     each_gene_info = gff_tomato[gff_tomato.iloc[:,8].str.contains(gene_id)]
#     each_gene_info["gene_id"] = gene_id
#     NLR_gene_info = pd.concat([NLR_gene_info, each_gene_info])
#
# NLR_gene_info = pd.DataFrame(NLR_gene_info)
# NLR_gene_info.to_csv("NLR_gene_info.csv")
# # print(NLR_gene_info.head())
# # print(NLR_gene_info.iloc[0,:])

# extract gene_distance info
NLR_gene_info = pd.read_csv("NLR_gene_info.csv", index_col=0)
position_data = NLR_gene_info.iloc[:, [0,3,4,9]]
position_data.columns = ["chr", "start", "end", "id"]

distances = []
def calc_distance(data_from_chr):
    tmp_position_data = data_from_chr.iloc[:, 1]
    tmp_position_data = np.sort(tmp_position_data)
    distance = np.diff(tmp_position_data, n= 1)
    distance = np.array(distance)
    distance = distance[distance < 50000]
    plt.hist(distance)
    plt.show()
    distances.extend(distance)

position_data.groupby(["chr"]).apply(calc_distance)
# calc_distance(position_data[position_data["chr"] == "SL3.0ch01"])

# distances = np.array(distances)
# print(distances)
# plt.hist(distances[distances < 50000])
# plt.hist(distances)
# plt.show()
