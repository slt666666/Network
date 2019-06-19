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

     Direction = dict(
        colorscale=[
            [0, "rgb(220, 220, 220)"],
            [0.2, "rgb(220, 220, 220)"],
            [0.2, "red"],
            [0.4, "red"],
            [0.4, "blue"],
            [0.6, "blue"],
            [0.6, "orange"],
            [0.8, "orange"],
            [0.8, "magenta"],
            [1.0, "magenta"]
        ],
        colorbar = dict(
            title = 'Gene direction',
            titleside = 'top',
            tickmode = 'array',
            tickvals = [0.4,1.2,2.0,2.8,3.6],
            tickwidth = 2,
            ticktext = ["nan", "TH ← ←", "HT → →", "HH → ←", "TT ← →"],
            ticks = 'outside',
            thickness = 20,
            len = 0.5
        )
     )

     ID = dict(
        colorscale="Reds",
        reversescale=False,
        # colorbar = dict(
        #     title = 'HMM_score',
        #     titleside = 'top',
        #     tickmode = 'array',
        #     tickvals = [0,5,10,15,20,25,30,35],
        #     ticktext = ["0","5","10","15","20","25","30","35"],
        #     ticks = 'outside'
        # )
     )

### for test set
class TestParams(Enum):

    clade_file = "sample_data/tomato_clade.csv"
    gff_file = "sample_data/tomato_NLR_gff.csv"
    expression_file = "sample_data/tomato_root_leaf_TPM.csv"
    coexpression_file = "sample_data/tomato_part_TPM.csv"
    MADA_file = "sample_data/tomato_MADA.csv"
    ID_file = "sample_data/tomato_ID.csv"
    tree_file = "sample_data/tomato_tree.nex"
    ordered_file = "sample_data/tomato_NLR_order.csv"
    distance_threshold = 30000
    phylogenetic_distance_threshold = 1

    Distance = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=None,
        threshold=distance_threshold,
        matrix_type="Distance",
        add_info=None
    )

    LogFC2 = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=None,
        threshold=distance_threshold,
        matrix_type="LogFC2",
        add_info=expression_file
    )

    Coexpression = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=None,
        threshold=distance_threshold,
        matrix_type="Coexpression",
        add_info=coexpression_file
    )

    MADA = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=None,
        threshold=distance_threshold,
        matrix_type="MADA",
        add_info=MADA_file
    )

    Direction = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=None,
        threshold=distance_threshold,
        matrix_type="Direction",
        add_info=None
    )

    ID = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=None,
        threshold=distance_threshold,
        matrix_type="ID",
        add_info=ID_file
    )

    PhyloDistance = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=ordered_file,
        nexus_tree=tree_file,
        threshold=phylogenetic_distance_threshold,
        matrix_type="Distance",
        add_info=None
    )

    PhyloLogFC2 = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=ordered_file,
        nexus_tree=tree_file,
        threshold=phylogenetic_distance_threshold,
        matrix_type="LogFC2",
        add_info=expression_file
    )

    PhyloCoexpression = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=ordered_file,
        nexus_tree=tree_file,
        threshold=phylogenetic_distance_threshold,
        matrix_type="Coexpression",
        add_info=coexpression_file
    )

    PhyloMADA = dict(
        clade_csv=clade_file,
        gff_csv=gff_file,
        ordered_csv=ordered_file,
        nexus_tree=tree_file,
        threshold=phylogenetic_distance_threshold,
        matrix_type="MADA",
        add_info=MADA_file
    )
