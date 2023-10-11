import pandas as pd
import csv

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Iterable, Union, Dict
from collections import OrderedDict
from itertools import starmap

from .utils import FilePath, ensure_path


SeqLike = Union[str, Seq, SeqRecord]

class Fasta(OrderedDict):

    def __init__(self, data: Iterable[SeqLike], id_prefix:str = '') -> None:
        self.id_prefix = id_prefix
        data = starmap(self.__value_check, enumerate(data))
        super().__init__(map(lambda x: (x.id, x), data))

    def __getitem__(self, __key: str) -> SeqRecord:
        return super().__getitem__(__key)

    def __setitem__(self, __key: str, __value: SeqRecord) -> None:
        return super().__setitem__(__key, __value)

    def __value_check(self, idx: int, seq: SeqLike) -> SeqRecord:
        if isinstance(seq, str):
            seq = Seq(seq)
        if isinstance(seq, Seq):
            seq = SeqRecord(
                seq,
                id=f'{idx}',
                description=''
            )
        if isinstance(seq, SeqRecord):
            seq.id = f'{self.id_prefix}{seq.id}'
            return seq
        raise ValueError('Element should be str, Seq or SeqRecord')

    def to_dict(self) -> Iterable[Dict]:
        """
        Iterate fasta sequence record.
        """
        for rid, record in self.items():
            yield {
                'id': rid,
                'sequence': str(record.seq),
                'description': record.description}
    
    def to_dataframe(self, id_as_index: bool = False) -> pd.DataFrame:
        """
        Read fasta sequence as `pd.DataFrame`.
        """
        df = pd.DataFrame(self.to_dict())
        if id_as_index:
            df.set_index('id', inplace=True)
        return df
    
    def to_fasta(self, path: FilePath, mkdir: bool = False):
        """
        Save as FASTA format.
        """
        save_fasta(self.values(), path=path, mkdir=mkdir)

    def to_csv(self, path: FilePath, mkdir: bool = False):
        """
        Save as CSV format.
        """
        path = ensure_path(path)
        if mkdir:
            path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, 'w') as f:
            writer = csv.DictWriter(
                f, 
                fieldnames=['id', 'sequence', 'description'],
                lineterminator='\n')
            writer.writeheader()
            writer.writerows(self.to_dict())


def read_fasta(path: FilePath) -> Fasta:
    if not isinstance(path, Path):
        path = Path(path)
    path = path.resolve().expanduser().absolute()
    return Fasta(SeqIO.parse(path, 'fasta'))


def save_fasta(sequences: Iterable[SeqRecord], path: FilePath, mkdir: bool = False):
    path = ensure_path(path)
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        SeqIO.write(sequences, f, 'fasta')


def df2fasta(df:pd.DataFrame,
             fasta_path: FilePath, 
             *seq_cols: str, 
             id_col: str = None, 
             mode: str = 'seperate', 
             sep: str = ''):
    """
    Convert a dataframe to a fasta file.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to be converted.
    fasta_path : str or Path
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
    if len(seq_cols) == 0:
        raise ValueError('seq_cols should not be empty')
    elif len(seq_cols) == 1:
        mode = 'joint'

    def _iter_seq():
        for index, row in df.iterrows():
            item_id = index if id_col is None else row[id_col]
            if mode == 'seperate':
                for seq_col in seq_cols:
                    seq_id = f"{item_id}_{seq_col}"
                    yield SeqRecord(Seq(row[seq_col]), id=seq_id, description='')
            elif mode == 'joint':
                seq = sep.join([row[seq_col] for seq_col in seq_cols])
                seq_id = str(item_id)
                yield SeqRecord(Seq(seq), id=seq_id, description='')
            else:
                raise ValueError(f"mode {mode} not supported")
    
    save_fasta(_iter_seq(), fasta_path)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
    csv2fasta_parser = subparsers.add_parser('csv2fasta')
    csv2fasta_parser.add_argument(
        '-i', '--input',
        required=True,
        type=Path,
        help='input csv file')
    csv2fasta_parser.add_argument(
        '-o', '--output',
        required=True,
        type=Path,
        help='output fasta file')
    csv2fasta_parser.add_argument(
        '-c', '--id_col',
        required=True,
        help='column name of the id')
    csv2fasta_parser.add_argument(
        '-s', '--seq_cols',
        required=True,
        nargs='+',
        help='column names of the sequences')
    csv2fasta_parser.add_argument(
        '--mode',
        default='seperate',
        choices=['seperate', 'joint'],
        help='mode of the conversion, \
        s means single sequence per entry, \
        j means joint sequence per entry')
    csv2fasta_parser.add_argument(
        '--sep',
        default='',
        help='separator of the joint sequence')
    
    fasta2csv_parser = subparsers.add_parser('fasta2csv')
    fasta2csv_parser.add_argument(
        '-i', '--input',
        required=True,
        type=Path,
        help='input fasta file')
    fasta2csv_parser.add_argument(
        '-o', '--output',
        required=True,
        type=Path,
        help='output csv file')
    args = parser.parse_args()

    if args.subcommand == 'csv2fasta':
        df = pd.read_csv(args.input)
        df2fasta(df, args.output, args.id_col, *args.seq_cols, args.mode, args.sep)

    elif args.subcommand == 'fasta2csv':
        fasta = read_fasta(args.input)
        fasta.to_csv(args.output)

    else:
        raise ValueError(f"subcommand {args.subcommand} not supported")
