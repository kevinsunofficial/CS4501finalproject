import os
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')


genomic = './sars_cov2_spike_aligned/sars_cov2_spike_aligned.fasta'
report = './sars_cov2_spike.csv'

sequences = {}
for seq_record in SeqIO.parse(genomic, 'fasta'):
    if len(seq_record.seq) == 3822:
        sequences[seq_record.id] = seq_record.seq

df = pd.read_csv(report)
df = df[df.Accession.isin(sequences.keys())]
df = df.sort_values(by='Accession')

accession_list = df.Accession.tolist()
len_acc = len(accession_list)

stat = np.zeros((len_acc,5,3822))
for i in tqdm(range(len_acc)):
    acc = accession_list[i]
    seq = sequences[acc]
    for j in range(3822):
        if seq[j] == 'A': 
            try: stat[i,0,j] = stat[i-1,0,j]+1
            except: stat[i,0,j] = 1
        elif seq[j] == 'T': 
            try: stat[i,1,j] = stat[i-1,1,j]+1
            except: stat[i,1,j] = 1
        elif seq[j] == 'C': 
            try: stat[i,2,j] = stat[i-1,2,j]+1
            except: stat[i,2,j] = 1
        elif seq[j] == 'G': 
            try: stat[i,3,j] = stat[i-1,3,j]+1
            except: stat[i,3,j] = 1
        elif seq[j] == '-': 
            try: stat[i,4,j] = stat[i-1,4,j]+1
            except: stat[i,4,j] = 1

fig, ax = plt.subplots()
ax.set_ylim(0,1.1)
palette = list((sns.color_palette("hls",5).as_hex()))
plt.style.use("seaborn")

x = [str(n) for n in range(500)]
width = 0.5

def animate(k):
    sums = np.sum(stat[k],axis=0)
    yA, yT, yC, yG, yn = stat[k,0], stat[k,1], stat[k,2], stat[k,3], stat[k,4]
    yA, yT, yC, yG, yn = yA/sums, yT/sums, yC/sums, yG/sums, yn/sums
    ax.bar(x, yA[-500:], label='A', width=width, color=palette[0])
    ax.bar(x, yT[-500:], label='T', width=width, color=palette[1], bottom=yA[-500:])
    ax.bar(x, yC[-500:], label='C', width=width, color=palette[2], bottom=yA[-500:]+yT[-500:])
    ax.bar(x, yG[-500:], label='G', width=width, color=palette[3], bottom=yA[-500:]+yT[-500:]+yC[-500:])
    ax.bar(x, yn[-500:], label='-', width=width, color=palette[4], bottom=yA[-500:]+yT[-500:]+yC[-500:]+yG[-500:])
    title = 'Distribution {}'.format(k+1) + ' up to date {}'.format(df[df.Accession==accession_list[k]].ReleaseDate.item())
    ax.set_title(title)
#     ax.legend(bbox_to_anchor=(1.04,0.5),loc='center left')
    ax.get_xaxis().set_visible(False)

anim = FuncAnimation(fig, animate, frames=25, interval=1, repeat=False)

anim.save('./test.gif', writer='imagemagick', fps=60)