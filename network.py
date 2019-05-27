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

# plotly.tools.set_credentials_file(username='slt666666', api_key='KEROY4DT9JraomPybSjQ')

# extract only "gene" list from gff
# gff_coffea = pd.read_table("Coffee/original_data/coffea_canephora.gff3", header=None, skiprows=1)
# gff_coffea = gff_coffea[gff_coffea.iloc[:, 2] == "gene"]
# # print(gff_tomato.head())
#
# # extract NLR gene_ids from data_AKI
# NLR_list = pd.read_csv("Coffee/original_data/coffee_NLR_id_class.csv")
# print(NLR_list)
#
# # extract gene info
# NLR_gene_info = pd.DataFrame()
# for gene_id in NLR_list["SequenceName"]:
#     each_gene_info = gff_coffea[gff_coffea.iloc[:, 8].str.contains(gene_id)]
#     each_gene_info["gene_id"] = gene_id
#     each_gene_info["class"] = NLR_list.loc[NLR_list["SequenceName"] == gene_id, "class"].values
#     NLR_gene_info = pd.concat([NLR_gene_info, each_gene_info])
#
# NLR_gene_info = pd.DataFrame(NLR_gene_info)
# NLR_gene_info.to_csv("Coffee/original_data/NLR_gene_info.csv")
# print(NLR_gene_info.head())
# print(NLR_gene_info.iloc[0,:])

# extract gene_distance info
NLR_gene_info = pd.read_csv("Coffee/original_data/NLR_gene_info.csv", index_col=0)
position_data = NLR_gene_info.iloc[:, [0,3,4,9]]
position_data.columns = ["chr", "start", "end", "id"]

distances = []
def calc_distance(data_from_chr):
    tmp_position_data = data_from_chr.iloc[:, 1]
    tmp_position_data = np.sort(tmp_position_data)
    distance = np.diff(tmp_position_data, n= 1)
    distance = np.array(distance)
    # distance = distance[distance < 50000]
    # plt.hist(distance)
    # plt.show()
    distances.extend(distance)

position_data.groupby(["chr"]).apply(calc_distance)
# calc_distance(position_data[position_data["chr"] == "SL3.0ch01"])

distances = np.array(distances)
# print(distances)
plt.hist(distances[distances < 100000])
# plt.hist(distances)
plt.show()
