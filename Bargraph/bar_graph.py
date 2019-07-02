import dendropy
import numpy as np
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import iplot, plot

import make_data


class BarGraph:

    def __init__(self, clade_info, bar_info, color_info, bar_type, color_type, main_plant=None, other_plant=None):
        self.clade_info = clade_info
        self.bar_info = bar_info
        self.color_info = color_info
        self.bar_type = bar_type
        self.color_type = color_type
        self.main_plant = main_plant
        self.other_plant = other_plant

    ### treat dataset
    def make_data(self):

        data_make = make_data.DataMake(self.clade_info, self.bar_info, self.color_info, self.bar_type, self.color_type, self.main_plant, self.other_plant)
        return data_make.make()

    ### draw graph
    def draw_bargraph(self):

        ### make data for plotly
        bar_data, bar_text, colors, annotations, title = self.make_data()

        ### bar graph datasetting
        trace = go.Bar(
            x=bar_data.index,
            y=bar_data,
            text=bar_text,
            marker=dict(
                color=colors
            ),
        )
        data = [trace]

        ### layout setting
        layout = go.Layout(
            margin=dict(
                l=130,
                b=160,
            ),
            width=1000,
            height=600,
            title=title,
            annotations=annotations
        )

        ### draw
        fig = go.Figure(
            data=data,
            layout=layout
        )
        plot(fig, filename='/Users/toshiyuki/Desktop/研究用/Network_data/Graphs/{}_{}.html'.format(self.main_plant, self.other_plant))


if __name__ == "__main__":
    for other_plant in ["Arabidopsis","Coffea","Kiwifruit","Lettuce","Pepper","Spinach","Sweetpotato","Tomato"]:
        plant = "Tomato"
        # other_plant = "Sweetpotato"
        clade_info = "/Users/toshiyuki/Desktop/研究用/Network_data/{}/original_data/{}_clade.csv".format(plant, plant)
        # bar_info = "/Users/toshiyuki/Desktop/研究用/Network_data/{}/Tomato{}.nwk".format(plant, plant)
        bar_info = "/Users/toshiyuki/Desktop/研究用/Network_data/Chrs/Chrs.nwk"
        color_info = "/Users/toshiyuki/Desktop/研究用/Network_data/{}/original_data/{}_MADA.csv".format(plant, plant)
        bar_graph = BarGraph(clade_info, bar_info, color_info, "Conserved", "MADA", plant, other_plant)
        # bar_graph = BarGraph("sample_data/tomato_clade.csv", "original_data/7species.tree", "sample_data/tomato_MADA.csv", "Conserved", "MADA", "Sweetpotato")
        # bar_graph = BarGraph("sample_data/tomato_clade.csv", "original_data/asterids.nwk", "sample_data/tomato_root_leaf_TPM.csv", "Conserved", "LogFC2", "Coffea")
        # bar_graph = BarGraph("sample_data/tomato_clade.csv", "original_data/tomato_root_leaf_TPM.csv", "sample_data/tomato_MADA.csv", "LogFC2", "MADA", "Arabidopsis")
        # bar_graph = BarGraph("sample_data/tomato_clade.csv", "original_data/tomato_MADA.csv", "sample_data/tomato_root_leaf_TPM.csv", "MADA", "LogFC2", "Arabidopsis")
        bar_graph.draw_bargraph()
