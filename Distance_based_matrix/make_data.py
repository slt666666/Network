import copy
import pandas as pd
import numpy as np
from scipy.spatial import distance


class DataMake:

    def __init__(self, clade_csv, gff_csv, threshold):
        self.clade_csv = clade_csv
        self.gff_csv = gff_csv
        self.threshold = threshold

    def make(self):

        ### clade information
        clade_info = pd.read_csv(self.clade_csv)
        clade_info.columns = ["id", "clade"]

        ### distance info
        gff_info = pd.read_csv(self.gff_csv)
        gff_info.columns = ["index", "chr", "source", "feature", "start", "end", "score", "strand", "frame", "group", "id"]

        ### merge data & add new ids
        position_data = pd.merge(gff_info, clade_info, on="id")
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
        z_data = copy.copy(dist)
        z_data[z_data==0] = self.threshold
        z_data[np.isnan(z_data)] = self.threshold
        z_data[z_data>self.threshold] = self.threshold

        ### make text of cells
        hovertext = list()
        for yi, yy in enumerate(ids[::-1]):
            hovertext.append(list())
            for xi, xx in enumerate(ids):
                hovertext[-1].append('1: {} {}<br />2: {} {}<br />Distance: {}'.format(yy, clades[::-1][yi], xx, clades[xi], dist[yi][xi]))

        return z_data, hovertext, position_data
