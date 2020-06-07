#%%
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import torch
import humanize
import math
np.random.seed(7) # fix random seed for reproducibility

from Bio import SeqIO,pairwise2
import Bio.SubsMat.MatrixInfo as matrices
from Bio.SubsMat import MatrixInfo as matlist

# %%
nucleotides = list(SeqIO.parse("Sars-Cov2-proteins.fasta", "fasta"))
proteins = list(SeqIO.parse("Sars-Cov2-nucleotides.fasta", "fasta"))
len(nucleotides),len(proteins)

#%%
nucleotide_len = [len(str(l.seq)) for l in nucleotides]
protein_len = [len(str(l.seq)) for l in proteins]


def plothist(data, bin_count):
    bins = np.linspace(math.ceil(min(data)),math.floor(max(data)),bin_count)
    plt.xlim([min(data)-5, max(data)+5])
    plt.hist(data, bins=bins, alpha=0.5)
    plt.ylabel('count')
    plt.show()

#%%
plothist(nucleotide_len, 50)
plothist(protein_len, 50)

#%%
def is_complete(seq):
    vals = seq.description.split('|')
    vals = [item.strip().lower() for item in vals]
    return (vals[5] == 'complete')
    
count = 0
for seq in nucleotides:
    if is_complete(seq):
        count+=1

print(count)

# %%
# Allocate memory
num_nu_seqs = len(nucleotides)
matrix = matlist.blosum62
num_clusters = 2

output = torch.zeros([num_nu_seqs,num_nu_seqs], dtype=torch.float32)
print("Allocated structure:",output.shape)
print("Allocated size:",humanize.naturalsize(sys.getsizeof(output)))


# %%
total_ops = num_nu_seqs**2
pbar = tqdm(total=total_ops)

for i in range(0,num_nu_seqs):
    for j in range(0,num_nu_seqs):
        a = pairwise2.align.globaldx(str(nucleotides[i].seq),str(nucleotides[j].seq), matrix)
        output[i][j] = a[0][2]
        pbar.update(1)

pbar.close()

# %%
