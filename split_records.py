#%%
from Bio import SeqIO

nucleotides = list(SeqIO.parse("Sars-Cov2-nucleotides.fasta", "fasta"))
proteins = list(SeqIO.parse("Sars-Cov2-proteins.fasta", "fasta"))

#%%
path = "nucleotides/"

#%%
print("Writing nucleotides...\n")
from tqdm import tqdm
for record in tqdm(nucleotides):
    handle = path + str(record.id) + ".fsa"
    SeqIO.write(record, handle,"fasta")

# %%
path = "proteins/"
print("Writing proteins...\n")
from tqdm import tqdm
for record in tqdm(proteins):
    handle = path + str(record.id) + ".fsa"
    SeqIO.write(record, handle,"fasta")