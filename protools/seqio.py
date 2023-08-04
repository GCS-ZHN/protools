import pandas as pd

from fasta import FASTA
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Iterable


def save_fasta(seq_data, fasta_file):
    """
    Save a fasta file.

    Parameters
    ----------
    seq_data : dict or iterable
        The sequence data to be saved.
    fasta_file : str
        The path of the fasta file to be saved.
    """
    if isinstance(seq_data, dict):
        save_fasta(seq_data.items(), fasta_file)
    
    elif isinstance(seq_data, Iterable):
        fasta_file = Path(fasta_file).resolve()
        if not fasta_file.parent.exists():
            fasta_file.parent.mkdir(parents=True)
        with FASTA(str(fasta_file)) as f:
            for seq_id, seq in seq_data:
                seq_record = SeqRecord(Seq(seq), id=seq_id, description='')
                f.add_seq(seq_record)
    else:
        raise TypeError('seq_data should be dict or iterable')


def df2fasta(df:pd.DataFrame,
             fasta_file:str, 
             id_col: str, 
             seq_cols: list, 
             mode: str = 'seperate', 
             sep: str = ''):
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
        The mode of the conversion. 'seperate' 
        for single sequence per entry,
        'joint' for joint sequence per entry.
        The default is 'seperate'.
    sep : str, optional
        The separator of the joint sequence. The default is ''.
        Only works when mode is 'joint'.
    """
    def _iter_seq():
        for _, row in df.iterrows():
            if mode == 'seperate':
                for seq_col in seq_cols:
                    seq_id = f"{row[id_col]}_{seq_col}"
                    yield seq_id, row[seq_col]
            elif mode == 'joint':
                seq = sep.join([row[seq_col] for seq_col in seq_cols])
                seq_id = row[id_col]
                yield seq_id, seq
            else:
                raise ValueError(f"mode {mode} not supported")
    
    save_fasta(_iter_seq(), fasta_file)



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='input csv file')
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='output fasta file')
    parser.add_argument(
        '-c', '--id_col',
        required=True,
        help='column name of the id')
    parser.add_argument(
        '-s', '--seq_cols',
        required=True,
        nargs='+',
        help='column names of the sequences')
    parser.add_argument(
        '--mode',
        default='seperate',
        choices=['seperate', 'joint'],
        help='mode of the conversion, \
        s means single sequence per entry, \
        j means joint sequence per entry')
    parser.add_argument(
        '--sep',
        default='',
        help='separator of the joint sequence')
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df2fasta(df, args.output, args.id_col, args.seq_cols, args.mode, args.sep)
