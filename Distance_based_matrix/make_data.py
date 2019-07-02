import copy
import dendropy
import pandas as pd
import numpy as np
from scipy.spatial import distance

### for genetic linkage distance
class DataMake:

    def __init__(self, clade_csv, gff_csv, ordered_csv, threshold, matrix_type, add_info):
        self.clade_csv = clade_csv
        self.gff_csv = gff_csv
        self.ordered_csv = ordered_csv
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

        ### make distance matrix
        dist = self.make_dist(starts, chrs)

        ### setting z_data for heatmap
        z_data = self.make_z_data(dist, matrix_type, position_data)

        ### make text of cells
        hovertext = self.make_hover_text(ids, clades, dist, matrix_type, position_data)

        return z_data, hovertext, position_data

    ### make distance matrix
    def make_dist(self, starts, chrs):

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

        return dist

    ### add additional information
    def make_additional_data(self, position_data, type, add_info):

        add_info = pd.read_csv(add_info)

        ### extract expression values of Root and Leaf
        if type == "LogFC2":
            add_info = add_info.loc[:, ["gene_short_na", "Root_TPM", "Leaf_TPM"]]
            add_info.columns = ["id", "Root_TPM", "Leaf_TPM"]
            position_data = pd.merge(position_data, add_info, on='id', how='left')

        ### extract expression values of each points
        elif type == "Coexpression":
            position_data = pd.merge(position_data, add_info, on='id', how='left')
            position_data = position_data.fillna(0)

        ### extract MADA HMM score
        elif type == "MADA":
            add_info = add_info.loc[:, ["ID", "HMM score", "Seq-F"]]
            add_info.columns = ["id", "HMM_score", "Seq-F"]
            position_data = pd.merge(position_data, add_info, on='id', how='left')

        elif type == "ID":
            position_data["Integrated"] = position_data.loc[:, "id"].isin(add_info["id"].values)

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

        ### matrix based on Coexpression values
        elif type == "Coexpression":
            co_expression_matrix = np.corrcoef(position_data.iloc[:, 12:22])
            self.origin_co_expression_matrix = copy.copy(co_expression_matrix)
            co_expression_matrix = np.flip(co_expression_matrix, 0)
            co_expression_matrix[z_data >= threshold] = 0
            z_data = co_expression_matrix

        ### matrix based on MADA HMM score
        elif type == "MADA":
            MADA = position_data["HMM_score"]
            MADA = MADA[::-1]
            MADA = MADA.fillna(0)
            MADA = np.tile(MADA, (position_data.shape[0],1)).T
            MADA[z_data >= threshold] = np.nan
            z_data = MADA

        ### matrix based on direction btw 2 genes
        elif type == "Direction":
            direction = position_data.loc[:, ["id", "strand"]]
            direction = direction.replace({"+":1, "-":2})
            direction = np.matrix(direction.loc[:, "strand"].values)
            direction = np.dot(direction.T, direction+1)
            direction = np.triu(direction)
            direction = direction + direction.T - np.diag(np.diag(direction))
            direction = np.flip(direction, 0)
            direction[z_data >= threshold] = 0
            direction = np.where(direction == 6, 1, direction)
            z_data = direction

        ### matrix based on Integrated domain
        elif type == "ID":
            ID = np.zeros(position_data.shape[0])
            ID[position_data["Integrated"]] = 1
            ID = ID[::-1]
            ID = np.tile(ID, (position_data.shape[0],1)).T
            ID[z_data >= threshold] = np.nan
            z_data = ID

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

                ### Co expression data
                elif type == "Coexpression":
                    hovertext[-1].append('1: {} {}<br />2: {} {}<br />correlation coef: {}<br />Distance: {}'.format(
                            yy,
                            clades[::-1][yi],
                            xx,
                            clades[xi],
                            self.origin_co_expression_matrix[position_data.shape[0] - 1 - yi][xi],
                            dist[yi][xi],
                        )
                    )

                ### MADA HMM score & MADA start position
                elif type == "MADA":
                    hovertext[-1].append('1: {} {}<br />2: {} {}<br />1: HMM_score: {} Seq-F: {}<br />2: HMM_score: {} Seq-F: {}<br />Distance: {}'.format(
                            yy,
                            clades[::-1][yi],
                            xx,
                            clades[xi],
                            position_data["HMM_score"].values[::-1][yi],
                            position_data["Seq-F"].values[::-1][yi],
                            position_data["HMM_score"].values[xi],
                            position_data["Seq-F"].values[xi],
                            dist[yi][xi],
                        )
                    )

                ### direction of each genes
                elif type == "Direction":
                    if np.isnan(dist[yi][xi]):
                        first_direction = "nan"
                        second_direction = "nan"
                    else:
                        first_direction = "→" if position_data["strand"].values[::-1][yi] == "+" else "←"
                        second_direction = "→" if position_data["strand"].values[xi] == "+" else "←"

                    if xi > position_data.shape[0] - 1 - yi:
                        hovertext[-1].append('1: {} {}<br />2: {} {}<br />direction: {} {}<br />Distance: {}'.format(
                                yy,
                                clades[::-1][yi],
                                xx,
                                clades[xi],
                                first_direction,
                                second_direction,
                                dist[yi][xi],
                            )
                        )
                    else:
                        hovertext[-1].append('1: {} {}<br />2: {} {}<br />direction: {} {}<br />Distance: {}'.format(
                                yy,
                                clades[::-1][yi],
                                xx,
                                clades[xi],
                                second_direction,
                                first_direction,
                                dist[yi][xi],
                            )
                        )

                    ### Integrated domain
                elif type == "ID":
                    hovertext[-1].append('1: {} {}<br />2: {} {}<br />1: Integrated Domain: {}<br />2: Integrated Domain: {}<br />Distance: {}'.format(
                            yy,
                            clades[::-1][yi],
                            xx,
                            clades[xi],
                            "○" if position_data["Integrated"].values[::-1][yi] else "×",
                            "○" if position_data["Integrated"].values[xi] else "×",
                            dist[yi][xi],
                        )
                    )

        return np.array(hovertext)


### for phylogenetic distance
class PhylogeneticDataMake(DataMake):

    def __init__(self, clade_csv, gff_csv, ordered_csv, nexus_tree, threshold, matrix_type, add_info):
        super().__init__(clade_csv, gff_csv, ordered_csv, threshold, matrix_type, add_info)
        self.nexus_tree = nexus_tree

    ### make phylogenetic distance matrix
    def make(self):

        ### define
        matrix_type = self.matrix_type
        add_info = self.add_info

        ### clade information
        clade_info = pd.read_csv(self.clade_csv)
        clade_info.columns = ["id", "clade"]

        ### distance info
        tree = dendropy.Tree.get(path=self.nexus_tree, schema="nexus", preserve_underscores=True)

        ### merge data & check other information
        position_data = clade_info

        if add_info != None and matrix_type != "Distance":
            position_data = self.make_additional_data(position_data, matrix_type, add_info)

        ### sort and add new ids, ordered genes if ordered_csv exist
        if self.ordered_csv == None:
            position_data = position_data.sort_values(["clade", "id"])
        else:
            NLR_order = pd.read_csv(self.ordered_csv, header=None)
            ordered_list = []
            for i in NLR_order.iloc[:, 0].values:
                if i in position_data["id"].values:
                    ordered_list.append(i)
            position_data.index = position_data["id"].values
            position_data = position_data.loc[ordered_list, :]

        position_data["id_clade"] = position_data["id"].str.cat(position_data["clade"], sep="-")

        ### each gff information
        ids = position_data["id"].values
        clades = position_data["clade"].values

        ### make distance matrix
        dist = self.make_dist(tree, position_data)

        ### setting z_data for heatmap
        z_data = self.make_z_data(dist, matrix_type, position_data)

        ### make text of cells
        hovertext = self.make_hover_text(ids, clades, dist, matrix_type, position_data)

        return z_data, hovertext, position_data

    ### make phylogenetic distance matrix
    def make_dist(self, tree, position_data):

        pdm = tree.phylogenetic_distance_matrix()
        labels = []
        dist = []
        for taxon1 in tree.taxon_namespace:
            labels.append(taxon1.label)
            each_rows = []
            for taxon2 in tree.taxon_namespace:
                weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
                each_rows.append(weighted_patristic_distance)
            dist.append(each_rows)
        dist = pd.DataFrame(dist)
        dist.index = labels
        dist.columns = labels
        dist = dist.loc[position_data.id, position_data.id]
        dist = np.array(dist)
        dist = np.flip(dist, 0)

        return dist
