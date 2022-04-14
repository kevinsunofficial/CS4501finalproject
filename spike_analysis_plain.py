import os
import csv
import pandas as pd

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


if __name__=='__main__':
    fastafn = 'sars_cov2_s_genomic.fasta'
    sequences = {}
    for seq_record in tqdm(SeqIO.parse(fastafn,'fasta'), desc='loading fasta file'):
        sequences[seq_record.id] = seq_record.seq
    print('loaded', len(sequences), 'sequences')

    csvfn = 'sars_cov2_s_report.csv'
    df = pd.read_csv(csvfn)
    df.head()

    refseq = sequences['NC_045512.2']

    s_sequences = [
        SeqRecord(
            refseq,
            id='NC_045512',
            name='S',
            description='surface glycoprotein, refseq'
        )
    ]

    for seq_id in tqdm(df.Accession.tolist(), desc='aligning'):
        seq2 = sequences[seq_id]
        alignment = pairwise2.align.globalms(refseq,seq2,2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)
        seqB = alignment[0].seqB
        
        s_sequences.append(
            SeqRecord(
                seqB,
                id=seq_id,
                name='S',
                description='surface glycoprotein'
            )
        )

    SeqIO.write(s_sequences, 'sars_cov2_s_aligned_genomic.fasta', 'fasta')