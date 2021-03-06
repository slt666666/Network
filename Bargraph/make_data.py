import copy
import dendropy
import numpy as np
import pandas as pd


class DataMake:

    def __init__(self, clade_info, bar_info, color_info, bar_type, color_type, main_plant, other_plant):
        self.clade_info = clade_info
        self.bar_info = bar_info
        self.color_info = color_info
        self.bar_type = bar_type
        self.color_type = color_type
        self.main_plant = main_plant
        self.other_plant = other_plant


    def make(self):

        bar_data = self.make_bar_data()
        colors, colors_data = self.make_colors(bar_data.index)
        bar_text = self.make_bar_text(bar_data, colors_data)
        annotations = self.make_annotations(bar_data)
        title = self.make_title()

        return bar_data, bar_text, colors, annotations, title


    ### make each bar values
    def make_bar_data(self):

        ### calculate the closest phylogenetic gene distances
        if self.bar_type == "Conserved":

            ### extract phylogenetic distance data from tree
            tree = dendropy.Tree.get(path=self.bar_info, schema="newick", preserve_underscores=True)
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

            plant_header = dict(
                Arabidopsis="AT",
                Coffea="Cc",
                Spinach="Spo",
                Sugarbeet="Bv",
                Ashtree="FRAEX",
                Sweetpotato="itf",
                Potato="PGSC",
                MonkeyFlower="Migut",
                Tobbaco="Niben",
                Pepper="CA",
                Kiwifruit="Achn",
                Tomato="Solyc",
                Lettuce="Lsat"
            )

            ### extract species specific (tomato/other) ids
            main_ids = distances.index[distances.index.str.contains(plant_header[self.main_plant])].values
            other_ids = distances.index[distances.index.str.contains(plant_header[self.other_plant])].values

            bar_data = distances.loc[main_ids, other_ids].min(axis=1)
            bar_data = bar_data.sort_values()

        ### get MADA HMM scores of all NLRs
        elif self.bar_type == "MADA":

            MADA = pd.read_csv(self.bar_info, index_col=0)
            MADA = MADA.loc[:, ["HMM score", "Seq-F"]]
            MADA.columns = ["HMM_score", "Seq-F"]
            MADA = MADA.sort_values("HMM_score")
            bar_data = MADA["HMM_score"]

        elif self.bar_type == "LogFC2":

            LogFC2 = pd.read_csv(self.bar_info, index_col=0)
            LogFC2 = LogFC2.loc[:, ["Root_TPM", "Leaf_TPM"]]
            LogFC2.columns = ["Root_TPM", "Leaf_TPM"]

            ### calculate logFC2
            logFC2_values = np.log2(LogFC2["Root_TPM"]) - np.log2(LogFC2["Leaf_TPM"])

            logFC2_values[logFC2_values > 10] = 10
            logFC2_values[logFC2_values < -10] = -10
            logFC2_values[np.isnan(logFC2_values)] = 0

            LogFC2["LogFC2"] = logFC2_values
            LogFC2 = LogFC2.sort_values("LogFC2")
            bar_data = LogFC2["LogFC2"]

        return bar_data

    ### make color list of each bar
    def make_colors(self, ordered_ids):

        if self.color_type == "Conserved":
            pass

        ### to colore MADA genes
        elif self.color_type == "MADA":

            base = 'cyan'
            colors = pd.Series([base]*len(ordered_ids.values))

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

            base = 'rgb(120, 120, 120)'
            colors = pd.Series([base]*len(ordered_ids.values))

            LogFC2 = pd.read_csv(self.color_info)
            LogFC2 = LogFC2.loc[:, ["gene_short_na", "Root_TPM", "Leaf_TPM"]]
            LogFC2.columns = ["id", "Root_TPM", "Leaf_TPM"]

            ### calculate logFC2
            colors_data = np.log2(LogFC2["Root_TPM"]) - np.log2(LogFC2["Leaf_TPM"])
            logFC2_values = colors_data
            LogFC2["LogFC2"] = logFC2_values
            colors_data = LogFC2

            ### adjust for color changes
            logFC2_values[logFC2_values > 10] = 10
            logFC2_values[logFC2_values < -10] = -10
            logFC2_values[np.isnan(logFC2_values)] = 0
            LogFC2["for_color"] = logFC2_values

            color_change_ids = ordered_ids.isin(LogFC2[LogFC2.loc[:, "for_color"] > 5].id.values)
            colors[color_change_ids] = 'rgb(165,0,38)'
            color_change_ids = ordered_ids.isin(LogFC2[(LogFC2.loc[:, "for_color"] <= 5) & (LogFC2.loc[:, "for_color"] > 0)].id.values)
            colors[color_change_ids] = 'rgb(244,109,67)'
            color_change_ids = ordered_ids.isin(LogFC2[LogFC2.loc[:, "for_color"] == 0].id.values)
            colors[color_change_ids] = 'rgb(120, 120, 120)'
            color_change_ids = ordered_ids.isin(LogFC2[(LogFC2.loc[:, "for_color"] < 0) & (LogFC2.loc[:, "for_color"] > -5)].id.values)
            colors[color_change_ids] = 'rgb(116,173,209)'
            color_change_ids = ordered_ids.isin(LogFC2[LogFC2.loc[:, "for_color"] < -5].id.values)
            colors[color_change_ids] = 'rgb(49,54,149)'

        return colors, colors_data

    ### make text to add each bar
    def make_bar_text(self, bar_data, colors_data):

        bar_text = []
        clade_info = pd.read_csv(self.clade_info)
        clade_info.columns = ["id", "clade"]

        for id in bar_data.index.values:

            ### extract clade information
            try:
                clade = clade_info.loc[clade_info["id"] == id, "clade"].values[0]
            except:
                clade = "NA"

            ### add MADA information
            if self.color_type == "MADA":
                try:
                    HMM_score = colors_data.loc[colors_data["id"] == id, "HMM_score"].values[0]
                    bar_text.append("HMM_score: {} <br />clade: {}".format(HMM_score, clade))
                except:
                    bar_text.append("clade: {}".format(clade))

            ### add LogFC2 information
            elif self.color_type == "LogFC2":
                try:
                    LogFC2 = colors_data.loc[colors_data["id"] == id, "LogFC2"].values[0]
                    bar_text.append("LogFC2: {} <br />clade: {}".format('{:.4f}'.format(LogFC2), clade))
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
            try:
                arrow_info = dict(
                    x=np.where(bar_data.index.values == k)[0][0],
                    y=bar_data.loc[k] if bar_data.loc[k] >= 0 else 0,
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
            except:
                pass

        ### mean annotations
        mean_info = dict(
            x = np.searchsorted(bar_data.values, bar_data.mean()),
            y = 0,
            text = "Mean distance"
        )
        annotations.append(mean_info)

        # if self.color_type == "MADA":
        #     pass
        # elif self.color_type == "LogFC2":
        #
        #     ### tribe info
        #     mean_info = dict(
        #         x = 5,
        #         y = 50,
        #         text = "Root",
        #         arrowhead = 7,
        #         arrowcolor = 'rgb(165,0,38)',
        #         ax=22,
        #         ay=0,
        #     )
        #     annotations.append(mean_info)
        #     mean_info = dict(
        #         x = 10,
        #         y = 50,
        #         text = "Leaf",
        #         arrowhead = 7,
        #         arrowcolor = 'rgb(49,54,149)',
        #         ax=22,
        #         ay=0,
        #     )
        #     annotations.append(mean_info)

        return annotations

    ### make title
    def make_title(self):

        if self.bar_type == "Conserved":
            title = 'Conserved btw {}/{}'.format(self.main_plant, self.other_plant)

        elif self.bar_type == "MADA":
            title = 'HMM score'

        elif self.bar_type == "LogFC2":
            title = 'LogFC2 of Root/Leaf TPM'

        if self.color_type == "Conserved":
            pass

        elif self.color_type == "MADA":
            title = title + " & colored MADA"

        elif self.color_type == "LogFC2":
            title = title + " & colored LogFC2 of Root/Leaf TPM"

        return title
