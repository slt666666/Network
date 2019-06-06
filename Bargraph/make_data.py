import dendropy
import numpy as np
import pandas as pd


class DataMake:

    def __init__(self, bar_info, color_info, bar_type, color_type):
        self.bar_info = bar_info
        self.color_info = color_info
        self.bar_type = bar_type
        self.color_type = color_type


    def make(self):

        bar_data = self.make_bar_data()
        colors, colors_data = self.make_colors(bar_data.index)
        bar_text = self.make_bar_text(bar_data, colors_data)
        annotations = self.make_annotations(bar_data)

        return bar_data, bar_text, colors, annotations


    ### make each bar values
    def make_bar_data(self):

        ### calculate the closest phylogenetic gene distances
        if self.bar_type == "Conserved":

            ### extract phylogenetic distance data from tree
            tree = dendropy.Tree.get(path=self.bar_info, schema="nexus")
            pdm = tree.phylogenetic_distance_matrix()

            labels = []
            distances = []
            for taxon1 in tree.taxon_namespace:
                labels.append(taxon1.label)
                each_rows = []
                for taxon2 in tree.taxon_namespace:
                    weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
                    each_rows.append(weighted_patristic_distance)
                distances.append(each_rows)

            distances = pd.DataFrame(distances)
            distances.index = labels
            distances.columns = labels

            ### extract species specific (tomato/other) ids
            tomato_ids = distances.index[distances.index.str.contains("Solyc")].values
            other_ids = distances.index[distances.index.str.contains("itf")].values

            bar_data = distances.loc[tomato_ids, other_ids].min(axis=1)
            bar_data = bar_data.sort_values()

        ### get MADA HMM scores of all NLRs
        elif self.bar_type == "MADA":

            MADA = pd.read_csv(self.bar_info, index_col=0)
            MADA = MADA.loc[:, ["HMM score", "Seq-F"]]
            MADA.columns = ["HMM_score", "Seq-F"]
            MADA = MADA.sort_values("HMM_score")
            bar_data = MADA

        elif self.bar_type == "LogFC2":
            pass

        return bar_data

    ### make color list of each bar
    def make_colors(self, ordered_ids):

        base = 'cyan'
        colors = pd.Series([base]*len(ordered_ids.values))

        if self.color_type == "Conserved":
            pass

        ### to colore MADA genes
        elif self.color_type == "MADA":

            MADA = pd.read_csv(self.color_info)
            MADA = MADA.loc[:, ["ID", "HMM score", "Seq-F"]]
            MADA.columns = ["id", "HMM_score", "Seq-F"]

            ### MADA genes
            MADA_ids = MADA[MADA.loc[:, "HMM_score"] >= 10].id.values
            color_change_ids = ordered_ids.isin(MADA_ids)
            colors[color_change_ids] = 'magenta'

            ### MADA Like genes
            MADA_ids = MADA[MADA.loc[:, "HMM_score"] < 10].id.values
            color_change_ids = ordered_ids.isin(MADA_ids)
            colors[color_change_ids] = 'orange'

            colors_data = MADA

        ### to color Root/Leaf specific
        elif self.color_type == "LogFC2":
            pass

        return colors, colors_data

    ### make text to add each bar
    def make_bar_text(self, bar_data, colors_data):

        bar_text = []

        if self.bar_type == "Conserved":
            clade_info = pd.read_csv("sample_data/tomato_clade.csv")
            clade_info.columns = ["id", "clade"]
            for id in bar_data.index.values:
                try:
                    clade = clade_info.loc[clade_info["id"] == id, "clade"].values[0]
                except:
                    clade = "NA"
                try:
                    HMM_score = colors_data.loc[colors_data["id"] == id, "HMM_score"].values[0]
                    bar_text.append("HMM_score: {} <br />clade: {}".format(HMM_score, clade))
                except:
                    bar_text.append("clade: {}".format(clade))

        return bar_text

    ### make some annotations
    def make_annotations(self, bar_data):

        annotations = []

        ### NRC & AtZAR1 dataset
        CNLs = dict(
            Solyc02g084890="CNL-ND-AtZAR1",
            Solyc07g053010="CNL-ND-AtZAR1",
            Solyc07g053020="CNL-ND-AtZAR1",
            Solyc01g090430="NRC1",
            Solyc03g005660="NRCn",
            Solyc04g007030="NRC4c",
            Solyc04g007050="NRC5",
            Solyc04g007060="NRC4b",
            Solyc04g007070="NRC4a",
            Solyc04g008150="NRC6",
            Solyc04g015210="NRC7",
            Solyc05g009630="NRC3",
            Solyc10g008220="NRC0",
            Solyc10g047320="NRC2",
        )

        for k, v in CNLs.items():
            arrow_info = dict(
                x=np.where(bar_data.index.values == k)[0][0],
                y=bar_data.loc[k],
                xref='x',
                yref='y',
                text=v,
                showarrow=True,
                arrowhead=1,
                arrowwidth=0.8,
                ax=0,
                ay=-100,
                hovertext=v
            )
            annotations.append(arrow_info)

        ### mean annotations
        mean_info = dict(
            x = np.searchsorted(bar_data.values, bar_data.mean()),
            y = 0,
            text = "Mean distance"
        )
        annotations.append(mean_info)

        return annotations
