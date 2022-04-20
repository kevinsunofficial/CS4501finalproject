import os
import csv
import pandas as pd

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

import multiprocessing

def aligning(seq_info):
    seq_id, seq2 = seq_info
    print(seq_id)
    alignment = pairwise2.align.globalms(refseq,seq2,2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)
    seqB = alignment[0].seqB
    return (seq_id, seqB)

if __name__=='__main__':
    fastafn = 'sars_cov2_spike.fasta'
    sequences = {}
    for seq_record in tqdm(SeqIO.parse(fastafn,'fasta'), desc='loading fasta file'):
        sequences[seq_record.id] = seq_record.seq
    print('loaded', len(sequences), 'sequences')

    csvfn = 'sars_cov2_spike.csv'
    df = pd.read_csv(csvfn)
    df.head()

    refseq = sequences['NC_045512.2']
    accession_list = df.Accession.tolist()
    accession_list.remove('NC_045512.2')
    align_seq = [(seq_id,sequences[seq_id]) for seq_id in tqdm(accession_list, desc='prepare alignemnt')]

    s_sequences = [
        SeqRecord(
            refseq,
            id='NC_045512.2',
            name='S',
            description='surface glycoprotein, refseq'
        )
    ]

    workers = multiprocessing.Pool(8)
    results = workers.map_async(aligning, align_seq).get()

    for seq_info in tqdm(results, desc='parsing result'):
        seq_id, seqB = seq_info
        s_sequences.append(
            SeqRecord(
                Seq(seqB),
                id=seq_id,
                name='S',
                description='surface glycoprotein'
            )
        )

    SeqIO.write(s_sequences, 'sars_cov2_s_aligned_genomic_parallel.fasta', 'fasta')
