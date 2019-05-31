from enum import Enum


class HeatMapParams(Enum):

     Distance = dict(
        colorscale="Reds",
        reversescale=True
     )

     LogFC2 = dict(
        colorscale="RdBu",
        reversescale=False,
        colorbar = dict(
            title = 'Log2FC<br />Root/Leaf<br />TPM values<br />&nbsp;',
            titleside = 'top',
            tickmode = 'array',
            tickvals = [-10,-5,0,5,10],
            ticktext = ["-10~<br />Leaf","-5","0","5","Root<br />10~"],
            ticks = 'outside'
        )
     )

     Coexpression = dict(
        colorscale="RdBu",
        reversescale=False,
        zmax=1,
        zmin=-1,
        colorbar = dict(
            title = 'correlation<br />coefficients',
            titleside = 'top',
            tickmode = 'array',
            tickvals = [-1,-0.5,0.,0.5,1,],
            ticktext = ["-1","-0.5","0","0.5","1"],
            ticks = 'outside'
        )
     )

     MADA = dict(
        colorscale="Reds",
        reversescale=False,
        colorbar = dict(
            title = 'HMM_score',
            titleside = 'top',
            tickmode = 'array',
            tickvals = [0,5,10,15,20,25,30,35],
            ticktext = ["0","5","10","15","20","25","30","35"],
            ticks = 'outside'
        )
     )
