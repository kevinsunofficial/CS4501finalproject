import os
import csv
import pandas as pd

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

import multiprocessing

fastafn = 'sars_cov2_s_genomic.fasta'
sequences = {}
for seq_record in tqdm(SeqIO.parse(fastafn,'fasta'), desc='loading fasta file'):
    sequences[seq_record.id] = seq_record.seq
print('loaded', len(sequences), 'sequences')

csvfn = 'sars_cov2_s_report.csv'
df = pd.read_csv(csvfn)
df.head()

refseq = sequences['NC_045512.2']
align_seq = [(seq_id,sequences[seq_id]) for seq_id in df.Accession.tolist()]

def aligning(seq_info):
    seq_id, seq2 = seq_info
    alignment = pairwise2.align.globalms(refseq,seq2,2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)
    return (seq_id, alignment[0].seqB)

workers = multiprocessing.Pool(4)
results = workers.map_async(aligning, align_seq).get()