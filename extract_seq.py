import pandas as pd
import numpy as np
import sys
from Bio import SeqIO


# ### extract fasta of NLR ids
# coffea_NLRs = pd.read_csv("Coffee/original_data/NLR_gene_info.csv")
#
# with open("Coffee/original_data/coffea_NLR_seq.fasta", mode='w') as f:
#
#     for record in SeqIO.parse('Coffee/original_data/coffea_pep.faa', 'fasta'):
#         id_part = record.id
#         desc_part = record.description
#         seq = record.seq
#
#         if id_part[:14] in coffea_NLRs["gene_id"].tolist():
#             print(id_part[:14])
#             f.write(">"+id_part[:14])
#             f.write("\n")
#             f.write(str(record.seq))
#             f.write("\n")


### cut based on K of GKTT
with open("Coffee/original_data/coffee_alignment_BN_extractAtZAR1_cutK.fasta", mode='w') as f:

    for record in SeqIO.parse('Coffee/original_data/coffee_alignment_BN_extractAtZAR1.fasta', 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = record.seq
        if str(record.seq)[50] == "K":
            # print(id_part)
            # print(str(record.seq)[49:53])
            f.write(">"+id_part)
            f.write("\n")
            f.write(str(record.seq))
            f.write("\n")
