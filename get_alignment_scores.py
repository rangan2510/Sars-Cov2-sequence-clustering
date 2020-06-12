#%%
import os
import sys
import time
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
import torch
import humanize
np.random.seed(7) # fix random seed for reproducibility

from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices

# %%
def is_complete_genome(seq):
    desc = seq.description
    # print(desc)
    if desc.find('complete genome') != -1:
        return True
    else:
        return False

# %%
nucleotides = list(SeqIO.parse("Sars-Cov2-nucleotides.fasta", "fasta"))
print(len(nucleotides),"nucleotide sequences avaialable.")

# %%
count = 0
culled_nucleotides = []
for seq in nucleotides:
    if is_complete_genome(seq):
        culled_nucleotides.append(seq)
        count+=1

print(count,"complete/partially complete genomes selected.")

#%%
rand_sample_set = 100
culled_nucleotides = random.sample(culled_nucleotides,rand_sample_set)
print("Randomly selected",rand_sample_set,"samples." )
# %%
num_nu_seqs = len(culled_nucleotides)
num_clusters = 2

output = torch.zeros([num_nu_seqs,num_nu_seqs], dtype=torch.float32)
print("Allocated structure:",output.shape)
print("Allocated size:",humanize.naturalsize(sys.getsizeof(output)))

#%%
total_ops = num_nu_seqs**2
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

# %%
pbar = tqdm(total=total_ops)
for i in range(0,num_nu_seqs):
    for j in range(0,num_nu_seqs):
        a = aligner.score(str(culled_nucleotides[i].seq),str(culled_nucleotides[j].seq))
        output[i][j] = a
        pbar.update(1)
pbar.close()
torch.save(output, 'output.pt') 