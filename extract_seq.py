import pandas as pd
import numpy as np
import sys
from Bio import SeqIO


tomato_NLRs = pd.read_csv("original_data/NLR_gene_info.csv")

with open("original_data/tomato_NLR_seq.fasta", mode='w') as f:

    for record in SeqIO.parse('original_data/ITAG3.2_proteins.fasta', 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = record.seq

        if id_part[:14] in tomato_NLRs["gene_id"].tolist():
            print(id_part[:14])
            f.write(">"+id_part[:14])
            f.write("\n")
            f.write(str(record.seq))
            f.write("\n")
