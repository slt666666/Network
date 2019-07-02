from plotly.offline import iplot, plot
import make_data

class ArrowGraph:

    def __init__(self, gff_info, clade_info, MADA_info, ID_info, threshold):
        self.gff_info = gff_info
        self.clade_info = clade_info
        self.MADA_info = MADA_info
        self.ID_info = ID_info
        self.threshold = threshold

    def make_data(self):
        data_make = make_data.DataMake(self.gff_info, self.clade_info, self.MADA_info, self.ID_info, self.threshold)
        return data_make.make_graph_info()

    def draw_arrows(self):
        fig = self.make_data()
        plot(fig, filename='Gene_cluster_MADA_IntegratedDomain.html')

if __name__ == "__main__":
    plant = "Arabidopsis"
    gff_info = "/Users/toshiyuki/Desktop/研究用/Network_data_old/{}/original_data/{}_NLR_gff.csv".format(plant, plant)
    clade_info = "/Users/toshiyuki/Desktop/研究用/Network_data_old/{}/original_data/{}_clade.csv".format(plant, plant)
    MADA_info = "/Users/toshiyuki/Desktop/研究用/Network_data_old/{}/original_data/{}_MADA.csv".format(plant, plant)
    ID_info = "/Users/toshiyuki/Desktop/研究用/Network_data/all/original_data/NLR_ID.csv"
    threshold = 30000
    arrow_graph = ArrowGraph(gff_info, clade_info, MADA_info, ID_info, threshold)
    arrow_graph.draw_arrows()
