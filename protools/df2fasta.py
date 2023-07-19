#!/usr/bin/env python
import pandas as pd
import os

from fasta import FASTA
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

def df2fasta(df, fasta_file, id_col: str, seq_cols: list, mode: str = 's', sep: str = ''):
    """
    Convert a dataframe to a fasta file.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to be converted.
    fasta_file : str
        The path of the fasta file to be saved.
    id_col : str
        The column name of the id.
    seq_cols : list
        The column names of the sequences.
    mode : str, optional
        The mode of the conversion. 's' for single sequence,
        'j' for joint sequence. The default is 's'.
    sep : str, optional
        The separator of the joint sequence. The default is ''.
        Only works when mode is 'j'.
    """
    fasta_file = os.path.abspath(fasta_file)
    with FASTA(fasta_file) as f:
        for _, row in df.iterrows():
            if mode == 's':
                for seq_col in seq_cols:
                    seq_id = f"{row[id_col]}_{seq_col}"
                    seq = Seq(row[seq_col])
                    seq_record = SeqRecord(seq, id=seq_id, description='')
                    f.add_seq(seq_record)
            elif mode == 'j':
                seq = Seq(sep.join([row[seq_col] for seq_col in seq_cols]))
                seq_id = row[id_col]
                seq_record = SeqRecord(seq, id=seq_id, description='')
                f.add_seq(seq_record)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='input csv file')
    parser.add_argument('-o', '--output', required=True, help='output fasta file')
    parser.add_argument('-c', '--id_col', required=True,  help='column name of the id')
    parser.add_argument('-s', '--seq_cols', required=True, nargs='+', help='column names of the sequences')
    parser.add_argument('-m', '--mode', default='s', help='mode of the conversion')
    parser.add_argument('--sep', default='', help='separator of the joint sequence')
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df2fasta(df, args.output, args.id_col, args.seq_cols, args.mode, args.sep)
