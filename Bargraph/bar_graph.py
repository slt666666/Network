import dendropy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

tree = dendropy.Tree.get(path="sample_data/TomatoPotato.nex", schema="nexus")
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
distances = np.array(distances)
distances = pd.DataFrame(distances)
distances.index = labels
distances.columns = labels

sweetpotato_ids = distances.index[distances.index.str.contains("itf")].values
tomato_ids = distances.index[distances.index.str.contains("Solyc")].values

closest_distance = distances.loc[tomato_ids, sweetpotato_ids].min(axis=1)
closest_distance = closest_distance.sort_values()

base = 'c'
newcolor = 'm'
colors = [base]*closest_distance.shape[0]
colors[5:100] = [newcolor]*95

plt.bar(range(closest_distance.shape[0]), closest_distance, color=colors)
plt.show()
