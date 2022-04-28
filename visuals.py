import numpy as np
import pandas as pd
from Bio import SeqIO

from tqdm import tqdm
import matplotlib.pyplot as plt

genomic = './sars_cov2_spike_aligned/sars_cov2_spike_aligned.fasta'
report = './sars_cov2_spike/sars_cov2_spike.csv'
sequences = {}
for seq_record in tqdm(SeqIO.parse(genomic, 'fasta')):
    if len(seq_record.seq) == 3822 and (set(seq_record.seq)<=set(['A','T','C','G','-','N','n'])):
        sequences[seq_record.id] = seq_record.seq

df = pd.read_csv(report)
df = df[df.Accession.isin(sequences.keys())]
df = df.sort_values(by='ReleaseDate', ascending=True)

accession_list = df.Accession.tolist()
accession_list.remove('NC_045512.2')
len_acc = len(accession_list)
refseq = sequences['NC_045512.2']

vcmap = {
    'alpha': 'grey',
    'Beta': 'red',
    'Gamma': 'orange',
    'Delta': 'green',
    'Lambda': 'blue',
    'Omicron': 'brown',
}

plt.figure(dpi=150)
for i in tqdm(range(len_acc)):
    acc = accession_list[i]
    seq = sequences[acc]
    var = df[df.Accession==acc].PangoClass.item()
    c = vcmap[var]
    a = 0.5 if var=='alpha' else 1
    for j in range(3822):
        if(seq[j]!=refseq[j]):
            plt.scatter(j,i,color=c,label=var.capitalize(),cmap='Pastel1',s=0.01, alpha=a)

plt.title('Mutation distribution')
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left")
plt.show()