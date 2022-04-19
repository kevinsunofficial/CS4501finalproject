import os
import csv
import pandas as pd

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
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
    accession_list = df.Accession.tolist()[:3]
    accession_list.remove('NC_045512.2')

    s_sequences = [
        SeqRecord(
            refseq,
            id='NC_045512.2',
            name='S',
            description='surface glycoprotein, refseq'
        )
    ]

    for seq_id in tqdm(accession_list, desc='aligning'):
        seq2 = sequences[seq_id]
        alignment1 = pairwise2.align.globalms(refseq[:1000],seq2[:1000],2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)
        alignment2 = pairwise2.align.globalms(refseq[1000:2000],seq2[1000:2000],2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)
        alignment3 = pairwise2.align.globalms(refseq[2000:3000],seq2[2000:3000],2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)
        alignment4 = pairwise2.align.globalms(refseq[3000:],seq2[3000:],2,-1,-10,-0.5,one_alignment_only=True,penalize_end_gaps=False)

        seqB = alignment1[0].seqB+alignment2[0].seqB+alignment3[0].seqB+alignment4[0].seqB
        
        s_sequences.append(
            SeqRecord(
                Seq(seqB),
                id=seq_id,
                name='S',
                description='surface glycoprotein'
            )
        )

    SeqIO.write(s_sequences, 'sars_cov2_s_aligned_genomic.fasta', 'fasta')