import pandas as pd
import tempfile
import csv

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Iterable, Union, Dict, Tuple, Any
from collections import OrderedDict
from itertools import starmap, product

from .utils import FilePath, ensure_path


SeqLike = Union[str, Seq, SeqRecord, Tuple[str, str]]

class Fasta(OrderedDict):

    def __init__(self, data: Iterable[SeqLike] = None, id_prefix:str = '') -> None:
        self.id_prefix = id_prefix
        if data is None:
            data = []
        data = starmap(self.__value_check, enumerate(data))
        super().__init__(map(lambda x: (x.id, x), data))

    def __getitem__(self, __key: Union[str, Iterable[str]]) -> Union[SeqRecord, 'Fasta']:
        try:
            return super().__getitem__(__key)
        except TypeError as e:
            if isinstance(__key, Iterable):
                sub_fasta = Fasta()
                sub_fasta.update(map(lambda x: (x, self[x]), __key))
                return sub_fasta
            raise e

    def __setitem__(self, __key: str, __value: SeqRecord) -> None:
        return super().__setitem__(__key, __value)

    def __value_check(self, idx: Any, seq: SeqLike) -> SeqRecord:
        if isinstance(seq, tuple):
            if len(seq) != 2:
                raise ValueError('Tuple should have 2 elements')
            if not all(isinstance(x, str) for x in seq):
                raise ValueError('Tuple should have 2 str elements')
            idx = seq[0]
            seq = seq[1]
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

    def __add__(self, other: 'Fasta') -> 'Fasta':
        """
        Concatenate two fasta objects.
        """
        new = self.copy()
        new.update(other)
        return new


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
            if isinstance(item_id, str):
                item_id = item_id.replace(' ', '_')
                item_id = item_id.replace('/', '_')
                item_id = item_id.replace(',', '_')
            if mode == 'seperate':
                for seq_col in seq_cols:
                    seq_id = f"{item_id}_{seq_col}"
                    if row[seq_col] == '':
                        continue
                    yield SeqRecord(Seq(row[seq_col]), id=seq_id, description='')
            elif mode == 'joint':
                seq = sep.join([row[seq_col] for seq_col in seq_cols])
                seq_id = str(item_id)
                yield SeqRecord(Seq(seq), id=seq_id, description='')
            else:
                raise ValueError(f"mode {mode} not supported")
    
    save_fasta(_iter_seq(), fasta_path)


def temp_fasta(path: FilePath, id_prefix: str = ''):
    path = ensure_path(path)
    fasta = read_fasta(path)
    id_map = dict()
    def _iter_seq():
        for i, (key, record) in enumerate(fasta.items()):
            record.id = f'{id_prefix}{i}'
            record.description = ''
            record.name = record.id
            id_map[record.id] = key
            yield record
    
    temp_file = tempfile.NamedTemporaryFile(
        prefix='temp_',
        suffix='.fasta',
        dir=path.parent)
    save_fasta(_iter_seq(), temp_file.name)
    return temp_file, id_map


def create_complex_seq(
        seq_id: str,
        *seqs: str,
        seq_description: str = '',
        linker: str = ':') -> SeqRecord:
    complex_seq = linker.join(seqs)
    return SeqRecord(
        Seq(complex_seq),
        id = seq_id,
        description=seq_description
    )


def cross_create(
        seq_records1: Iterable[SeqRecord],
        seq_records2: Iterable[SeqRecord],
        linker: str = ':'
) -> Iterable[SeqRecord]:
    for seq1, seq2 in product(seq_records1, seq_records2):
        complex_id = f'{seq1.id}-{seq2.id}_complex'
        yield create_complex_seq(
            complex_id,
            str(seq1.seq),
            str(seq2.seq),
            linker=linker)


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

    complex_parser = subparsers.add_parser('complex')
    complex_parser.add_argument('--seqs1', '-i1', type=Path, required=True)
    complex_parser.add_argument('--seqs2', '-i2', type=Path, required=True)
    complex_parser.add_argument('--output', '-o', type=Path, required=True)
    complex_parser.add_argument('--linker', '-l', default='')

    args = parser.parse_args()

    if args.subcommand == 'csv2fasta':
        df = pd.read_csv(args.input)
        df2fasta(df, args.output, args.id_col, *args.seq_cols, args.mode, args.sep)

    elif args.subcommand == 'fasta2csv':
        fasta = read_fasta(args.input)
        fasta.to_csv(args.output)

    elif args.subcommand == 'complex':
        seqs1 = read_fasta(args.seqs1)
        seqs2 = read_fasta(args.seqs2)

        save_fasta(
            cross_create(seqs1.values(), seqs2.values(), args.linker),
            args.output)

    else:
        raise ValueError(f"subcommand {args.subcommand} not supported")
