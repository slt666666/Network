import random
import make_data
import pandas as pd
import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import iplot, plot
from scipy.spatial import distance


class DistanceBasedMatrix:

    def __init__(self, clade_csv, gff_csv, threshold):
        self.clade_csv = clade_csv
        self.gff_csv = gff_csv
        self.threshold = threshold

    def draw_heatmap(self, filename):

        ### make dataset
        data_make = make_data.DataMake(self.clade_csv, self.gff_csv, self.threshold)
        dist, hovertext, position_data = data_make.make()

        id_clades = position_data["id_clade"]
        gene_num = position_data.shape[0]
        clades = position_data["clade"]

        ### base of heatmap
        trace = go.Heatmap(
            x=id_clades,
            y=id_clades.values[::-1],
            z=dist,
            colorscale="Reds",
            reversescale=True,
            hoverinfo='text',
            text=hovertext
        )
        data = [trace]

        ### add separate lines of clade
        lines, colorlist = [], []
        start = -0.5

        for i, each_clade in enumerate(clades.unique()):

            ### color list of clade line
            colorlist.append("rgb"+str(tuple(np.random.random(size=3)*256)))

            start += clades.value_counts()[each_clade]

            ### vertical lines
            lines.append(
                dict(
                    type="line",
                    xref="x",
                    yref="y",
                    x0=start,
                    y0=-0.5,
                    x1=start,
                    y1=gene_num-0.5,
                    opacity=0.5,
                    line=dict(
                        color=colorlist[i],
                        width=0.5
                    )
                )
            )
            lines.append(
                dict(
                    type="line",
                    xref="x",
                    yref="y",
                    x0=-0.5,
                    y0=gene_num-1-start,
                    x1=gene_num-0.5,
                    y1=gene_num-1-start,
                    opacity=0.75,
                    line=dict(
                        color=colorlist[i],
                        width=0.5
                    )
                )
            )

        ### setting layout
        layout = go.Layout(
            margin=dict(l=130, t=160),
            width=750,
            height=750,
            autosize=False,
            xaxis=dict(
                    mirror="allticks",
                    side="top",
                    tickfont=dict(size=8)
            ),
            yaxis=dict(
                tickfont=dict(size=8)
            ),
            shapes=lines
        )

        fig = go.Figure(
            data=data,
            layout=layout
        )

        ### draw
        plot(fig, filename=filename)

if __name__ == "__main__":
    clade_file = "original_data/tomato_clade.csv"
    gff_file = "original_data/tomato_NLR_gff.csv"
    threshold = 30000

    Matrix = DistanceBasedMatrix(clade_file, gff_file, threshold)
    Matrix.draw_heatmap("sample.html")