import dendropy
import numpy as np
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import iplot, plot

import make_data


class BarGraph:

    def __init__(self, bar_info, color_info, bar_type, color_type):
        self.bar_info = bar_info
        self.color_info = color_info
        self.bar_type = bar_type
        self.color_type = color_type

    ### treat dataset
    def make_data(self):

        data_make = make_data.DataMake(self.bar_info, self.color_info, self.bar_type, self.color_type)
        return data_make.make()

    ### draw graph
    def draw_bargraph(self):

        ### make data for plotly
        bar_data, bar_text, colors, annotations = self.make_data()

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
            title='Conserved btw Tomato/Sweetpotato & colored MADA',
            annotations=annotations
        )

        ### draw
        fig = go.Figure(
            data=data,
            layout=layout
        )
        plot(fig, filename='sample.html')


if __name__ == "__main__":
    bar_graph = BarGraph("original_data/sweetpotato_tree.nex", "sample_data/tomato_MADA.csv", "Conserved", "MADA")
    bar_graph.draw_bargraph()
