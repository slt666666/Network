import numpy as np
import pandas as pd
import colorlover as cl

from scipy.spatial import distance

import plotly.plotly as py
import plotly.graph_objs as go
from plotly import tools


class DataMake:

    def __init__(self, gff_info, clade_info, MADA_info, ID_info, threshold):
        self.gff_info = gff_info
        self.clade_info = clade_info
        self.MADA_info = MADA_info
        self.ID_info = ID_info
        self.threshold = threshold

    def calc_distance(self, position_data):

        ### make distance matrix
        chrs = position_data["chr"].values
        dim = position_data.shape[0]
        start_position = position_data.loc[:, "start"]
        end_position = position_data.loc[:, "end"]
        start_matrix = np.repeat(np.array(start_position), dim).reshape(dim, dim).T
        end_matrix = np.repeat(np.array(end_position), dim).reshape(dim, dim)

        distance = start_matrix - end_matrix
        distance = np.triu(distance, 1)
        distance = distance + distance.T

        distance = pd.DataFrame(distance)
        distance.index = chrs
        distance.columns = chrs

        ### set np.nan to diffrent chromosome gene pair
        for each_chr in distance.columns:
            distance.loc[distance.index != each_chr, each_chr] = np.nan

        ### fix some misordered gene distance
        mistake_index = np.where(distance < 0)
        for i in range(len(mistake_index[0])):
            first = position_data.iloc[mistake_index[0][i], :]
            second = position_data.iloc[mistake_index[1][i], :]
            if first.start < second.start:
                distance.iloc[mistake_index[0][i], mistake_index[1][i]] = second.start - first.end
            else:
                distance.iloc[mistake_index[0][i], mistake_index[1][i]] = first.start - second.end

        distance = np.array(distance)
        return distance

    def make_position_data(self):

        threshold = self.threshold

        ### extract clade information
        clade_info = pd.read_csv(self.clade_info)
        clade_info.columns = ["id", "clade"]

        ### extract gene_distance info
        NLR_gene_info = pd.read_csv(self.gff_info, index_col=0)
        position_data = NLR_gene_info.iloc[:, [0,3,4,6,9]]
        position_data.columns = ["chr", "start", "end", "direction", "id"]

        ### MADA information
        MADA_info = pd.read_csv(self.MADA_info)
        MADA_info = MADA_info.loc[:, ["ID", "HMM score", "Seq-F"]]
        MADA_info.columns = ["id", "HMM_score", "Seq_F"]

        ### ID information
        ID_info = pd.read_csv(self.ID_info, header=None)
        ID_info = position_data["id"].isin(ID_info.iloc[:, 0])
        position_data["Integrated"] = ID_info

        position_data = pd.merge(position_data, clade_info, on='id')
        position_data = pd.merge(position_data, MADA_info, on='id', how='left')

        position_data = position_data.sort_values(["id", "start"])

        return position_data

    def make_color_list(self, position_data):

        ### make color list for each clade
        clade_list = position_data["clade"].unique()
        color_list = cl.scales['11']['div']['Spectral']
        color_list = cl.interp(color_list, len(clade_list))
        color_list = dict(zip(clade_list, color_list))

        return color_list

    def search_cluster(self, base, cluster_ids, rest_data, all_clusters):

        ### recursion methods
        ### search clustered gene with base gene
        threshold = self.threshold
        cluster_ids.append(base)

        ### end
        if rest_data.shape[0] <= 1:
            all_clusters.append(cluster_ids)

        ### start search
        else:

            ### if there are some clustered genes
            if rest_data.loc[base, :].min() < threshold:
                rm_ids = rest_data.loc[:, rest_data.loc[base, :] < threshold].columns.values
                cluster_ids.extend(rm_ids[:-1])
                rest_data = rest_data.drop(base)
                rest_data = rest_data.drop(base, axis=1)
                rest_data = rest_data.drop(rm_ids[:-1], axis=0)
                rest_data = rest_data.drop(rm_ids[:-1], axis=1)
                self.search_cluster(rm_ids[-1], cluster_ids, rest_data, all_clusters)

            ### if not, finish one cluster
            else:
                all_clusters.append(cluster_ids)
                rest_data = rest_data.drop(base, axis=0)
                rest_data = rest_data.drop(base, axis=1)
                self.search_cluster(rest_data.columns[0], [], rest_data, all_clusters)

        return all_clusters

    def make_clusters(self, position_data):

        threshold = self.threshold

        ### calculate distance btw genes
        dist = pd.DataFrame(self.calc_distance(position_data))
        dist.index = position_data["chr"].values
        dist.columns = position_data["chr"].values

        for i in range(dist.shape[0]):
            diff_chr = dist.columns != dist.index[i]
            dist.iloc[i, :].loc[diff_chr] = threshold

        dist = np.array(dist)
        dist[dist==0] = threshold
        distance_data = pd.DataFrame(dist)
        distance_data.index = position_data["id"].values
        distance_data.columns = position_data["id"].values

        ### make all clusters and all genes in cluster
        all_clusters = self.search_cluster(position_data["id"].values[0], [], distance_data, [])

        ### extract cluster with more than 2 genes
        clusters = []
        for cluster in all_clusters:
            if len(cluster) > 1:
                clusters.append(cluster)

        return clusters

    def make_graph_info(self):

        position_data = self.make_position_data()
        clusters = self.make_clusters(position_data)
        # color_list = self.make_color_list(position_data)

        ### decide xaxis range and sort clusters based on number of genes
        widest_dist = 0
        element_num = []
        for cluster in clusters:
            element_num.append(len(cluster))
            cluster_data = position_data.loc[position_data["id"].isin(cluster), :]
            base_line_start = cluster_data["start"].min()
            base_line_end = cluster_data["end"].max()
            if base_line_end - base_line_start > widest_dist:
                widest_dist =  base_line_end - base_line_start

        ### make subplots
        fig = tools.make_subplots(rows=len(clusters), cols=1)

        lines = []
        chrs = []
        for i, cluster_index in enumerate(np.argsort(element_num)[::-1]):
            cluster = clusters[cluster_index]
            cluster_data = position_data.loc[position_data["id"].isin(cluster), :]
            base_line_start = cluster_data["start"].min()
            chrs.append(cluster_data["chr"].unique()[0])

            ### make color_list
            color_list = ["red", "blue", "green", "orange", "purple"]
            color_list = dict(zip(cluster_data["clade"].unique(), color_list))

            ### string of axis number
            if i == 0:
                axis_num = ""
            else:
                axis_num = i+1

            ### draw base line of each cluster
            lines.append(
                dict(
                    type="line",
                    xref="x{}".format(axis_num),
                    yref="y{}".format(axis_num),
                    x0=base_line_start,
                    y0=i+100,
                    x1=base_line_start + widest_dist,
                    y1=i+100,
                    line=dict(
                        color="black"
                    ),
                    layer="below"
                )
            )

            ### draw each gene arrow and line and MADA and ...etc
            gene_arrow_x = []
            gene_symbol = []
            gene_texts = []
            gene_color = []
            MADA_x = []
            MADA_size = []
            MADA_symbol = []
            Integrated_x = []
            for j in range(cluster_data.shape[0]):

                ### check duplicate position (Not Yet)
                each_gene_data = cluster_data.iloc[j, :]

                lines.append(
                    dict(
                        type="line",
                        xref="x{}".format(axis_num),
                        yref="y{}".format(axis_num),
                        x0=each_gene_data.start,
                        y0=i+100,
                        x1=each_gene_data.end,
                        y1=i+100,
                        line=dict(
                            color=color_list[each_gene_data.clade],
                            width=15
                        ),
                        layer="below"
                    )
                )
                gene_color.append(color_list[each_gene_data.clade])
                gene_text = "{}<br>{} {}-{}<br>{}".format(each_gene_data.id, each_gene_data.chr, each_gene_data.start, each_gene_data.end, each_gene_data.clade)

                ### check direction of gene
                if each_gene_data.direction == "+":
                    gene_arrow_x.append(each_gene_data.end)
                    gene_symbol.append(8)
                elif each_gene_data.direction == "-":
                    gene_arrow_x.append(each_gene_data.start)
                    gene_symbol.append(7)

                ### check MADA
                if each_gene_data.HMM_score > 0:
                    MADA_x.append(each_gene_data.end if each_gene_data.direction == "+" else each_gene_data.start)
                    MADA_size.append(each_gene_data.HMM_score * 2 / 3)
                    MADA_symbol.append(8 if each_gene_data.direction == "+" else 7)
                    gene_text += "<br>HMM_score:{}, Seq-F:{}".format(each_gene_data.HMM_score, each_gene_data.Seq_F)

                gene_texts.append(gene_text)

                ### check integrated domain
                if each_gene_data.Integrated:
                    Integrated_x.append(each_gene_data.start if each_gene_data.direction == "+" else each_gene_data.end)

            ### draw triangle of arrows
            genes = go.Scatter(
                x=gene_arrow_x,
                y=np.ones(len(gene_arrow_x)) * (i+100),
                xaxis="x{}".format(axis_num),
                yaxis="y{}".format(axis_num),
                marker=dict(
                    color=gene_color,
                    symbol=gene_symbol,
                    size=30,
                ),
                showlegend=False,
                mode="markers",
                text=gene_texts,
                hoverinfo="x+text",
            )
            fig.append_trace(genes, i+1, 1)

            ### draw triangle of MADAs
            if len(MADA_x) > 0:
                MADAs = go.Scatter(
                    x=MADA_x,
                    y=np.ones(len(MADA_x)) * (i+100),
                    xaxis="x{}".format(axis_num),
                    yaxis="y{}".format(axis_num),
                    marker=dict(
                        color="black",
                        symbol=MADA_symbol,
                        size=MADA_size,
                    ),
                    showlegend=False,
                    mode="markers",
                    hoverinfo="none",
                )
                fig.append_trace(MADAs, i+1, 1)

            ### draw Integrated domains
            if len(Integrated_x) > 0:
                IDs = go.Scatter(
                    x=Integrated_x,
                    y=np.ones(len(Integrated_x)) * (i+100),
                    xaxis="x{}".format(axis_num),
                    yaxis="y{}".format(axis_num),
                    marker=dict(
                        color="black",
                        symbol=3,
                        size=15,
                    ),
                    showlegend=False,
                    mode="markers",
                    hoverinfo="none",
                )
                fig.append_trace(IDs, i+1, 1)

        fig['layout'].update(
            shapes=lines,
            height=len(clusters)*150
        )

        for i in range(len(clusters)):
            axis_num = i+1
            fig['layout']['xaxis{}'.format(axis_num)].update(
                showgrid=False
            )
            fig['layout']['yaxis{}'.format(axis_num)].update(
                title=dict(
                    text="{}".format(chrs[i])
                ),
                showgrid=False,
                showticklabels=False
            )

        return fig
