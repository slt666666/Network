import copy
import pandas as pd
import numpy as np
from scipy.spatial import distance


class DataMake:

    def __init__(self, clade_csv, gff_csv, threshold, matrix_type, add_info):
        self.clade_csv = clade_csv
        self.gff_csv = gff_csv
        self.threshold = threshold
        self.matrix_type = matrix_type
        self.add_info = add_info


    def make(self):

        ### define
        matrix_type = self.matrix_type
        add_info = self.add_info

        ### clade information
        clade_info = pd.read_csv(self.clade_csv)
        clade_info.columns = ["id", "clade"]

        ### distance info
        gff_info = pd.read_csv(self.gff_csv)
        gff_info.columns = ["index", "chr", "source", "feature", "start", "end", "score", "strand", "frame", "group", "id"]

        ### merge data & check other information
        position_data = pd.merge(gff_info, clade_info, on="id")

        if add_info != None and matrix_type != "Distance":
            position_data = self.make_additional_data(position_data, matrix_type, add_info)

        ### sort and add new ids
        position_data = position_data.sort_values(["clade", "id"])
        position_data["id_clade"] = position_data["id"].str.cat(position_data["clade"], sep="-")

        ### each gff information
        ids = position_data["id"].values
        chrs = position_data["chr"].values
        clades = position_data["clade"].values
        starts = position_data["start"].values

        ### calculate distance all combination
        dist = np.stack([np.array(starts), np.zeros(len(starts))], 1)
        dist = distance.cdist(dist, dist, metric="euclidean")
        dist = pd.DataFrame(dist)

        dist.index = chrs
        dist.columns = chrs

        ### set np.nan to diffrent chromosome gene pair
        for each_chr in dist.columns:
            dist.loc[dist.index != each_chr, each_chr] = np.nan
        dist = np.array(dist)
        dist = np.flip(dist, 0)

        ### setting z_data for heatmap
        z_data = self.make_z_data(dist, matrix_type, position_data)

        ### make text of cells
        hovertext = self.make_hover_text(ids, clades, dist, matrix_type, position_data)

        return z_data, hovertext, position_data


    ### add additional information
    def make_additional_data(self, position_data, type, add_info):

        add_info = pd.read_csv(add_info)

        ### extract expression values od Root and Leaf
        if type == "LogFC2":
            add_info = add_info.loc[:, ["gene_short_na", "Root_TPM", "Leaf_TPM"]]
            add_info.columns = ["id", "Root_TPM", "Leaf_TPM"]
            position_data = pd.merge(position_data, add_info, on='id')

        elif type == "Coexpression":
            pass
        elif type == "MADA":
            pass

        return position_data

    ### make z_data for plotly heatmap
    def make_z_data(self, dist, type, position_data):

        threshold = self.threshold

        ### matrix based on distance (default)
        z_data = copy.copy(dist)
        z_data[z_data==0] = threshold
        z_data[np.isnan(z_data)] = threshold
        z_data[z_data > threshold] = threshold

        ### matrix based on Root/Leaf logFC2 values
        if type == "LogFC2":
            logFC2 = np.log2(position_data["Root_TPM"]) - np.log2(position_data["Leaf_TPM"])
            logFC2 = logFC2.values[::-1]
            logFC2[logFC2 > 10] = 10
            logFC2[logFC2 < -10] = -10
            logFC2[np.isnan(logFC2)] = 0
            logFC2 = np.tile(logFC2, (position_data.shape[0],1)).T
            logFC2[z_data >= threshold] = 0
            z_data = logFC2

        else:
            pass

        return z_data

    ### make each cell text
    def make_hover_text(self, ids, clades, dist, type, position_data):

        hovertext = list()
        for yi, yy in enumerate(ids[::-1]):
            hovertext.append(list())
            for xi, xx in enumerate(ids):

                ### basic hover text
                if type == "Distance":
                    hovertext[-1].append('1: {} {}<br />2: {} {}<br />Distance: {}'.format(
                            yy,
                            clades[::-1][yi],
                            xx,
                            clades[xi],
                            dist[yi][xi]
                        )
                    )

                ### Root & Leaf expression data
                elif type == "LogFC2":
                    hovertext[-1].append('1: {} {}<br />2: {} {}<br />1: Root: {} Leaf: {}<br />2: Root: {} Leaf: {}<br />Distance: {}'.format(
                            yy,
                            clades[::-1][yi],
                            xx,
                            clades[xi],
                            '{:.4f}'.format(position_data["Root_TPM"].values[::-1][yi]),
                            '{:.4f}'.format(position_data["Leaf_TPM"].values[::-1][yi]),
                            '{:.4f}'.format(position_data["Root_TPM"].values[xi]),
                            '{:.4f}'.format(position_data["Leaf_TPM"].values[xi]),
                            dist[yi][xi],
                        )
                    )

                else:
                    pass



        return hovertext
