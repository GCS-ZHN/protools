import pandas as pd
import tempfile
import csv

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Iterable, Union, Dict, Tuple, Optional
from collections import OrderedDict
from itertools import product

from .utils import ensure_fileio, ensure_path
from .typedef import FilePathType, FilePathOrIOType, SeqLikeType


class Fasta(OrderedDict):
    """
    Sequence record dictionary
    and FASTA format processor.
    
    Parameters
    ----------
    data : Iterable[Tuple[str, SeqLike]], optional
        The initial data of the Fasta object.
        The default is None.
    **kwargs : Dict[str, SeqLike]
        The initial data of the Fasta object.
        The default is None.

    Examples
    ----------
    >>> from protools.seqio import Fasta
    >>> fasta = Fasta([('seq1', 'ATCG'), ('seq2', 'ATCG')])
    >>> fasta['seq1']
    SeqRecord(seq=Seq('ATCG'), id='seq1', name='seq1', description='', dbxrefs=[])
    """
    def __init__(self, 
                 data: Optional[Iterable[Tuple[str, SeqLikeType]]] = None, 
                 **kwargs: Dict[str, SeqLikeType]) -> None:
        if data is None:
            data = []
        super().__init__(data, **kwargs)

    def __getitem__(self, __key: Union[str, Iterable[str]]) -> Union[SeqRecord, 'Fasta']:
        try:
            return super().__getitem__(__key)
        except TypeError as e:
            if isinstance(__key, Iterable):
                sub_fasta = Fasta()
                sub_fasta.update(map(lambda x: (x, self[x]), __key))
                return sub_fasta
            raise e

    def __setitem__(self, __key: str, __value: SeqLikeType) -> None:
        __key = str(__key)
        if isinstance(__value, str):
            __value = Seq(__value)
        if isinstance(__value, Seq):
            __value = SeqRecord(
                __value,
                id=__key,
                description='',
                name=__key
            )
        if isinstance(__value, SeqRecord):
            assert __value.id == __key, f"Mismatch id and SeqRecord: {__key} != {__value.id}"
            if hasattr(self, '_binded_file'):
                save_fasta([__value], self._binded_file, **self._binded_kwargs)
            return super().__setitem__(__key, __value)
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
    
    def to_fasta(self, path: FilePathOrIOType, mkdir: bool = False, **kwargs):
        """
        Save as FASTA format.
        """
        save_fasta(self.values(), path=path, mkdir=mkdir, **kwargs)

    def to_csv(self, path: FilePathOrIOType, mkdir: bool = False):
        """
        Save as CSV format.
        """
        f, need_close = ensure_fileio(path, 'w')
        if mkdir:
            path.parent.mkdir(parents=True, exist_ok=True)
        writer = csv.DictWriter(
            f, 
            fieldnames=['id', 'sequence', 'description'],
            lineterminator='\n')
        writer.writeheader()
        writer.writerows(self.to_dict())
        if need_close:
            f.close()
        else:
            f.flush()

    def __add__(self, other: 'Fasta') -> 'Fasta':
        """
        Concatenate two fasta objects.
        """
        new = self.copy()
        new.update(other)
        return new
    
    def unique(self) -> 'Fasta':
        uniqued = Fasta()
        seqs = set()
        for rid, record in self.items():
            seq = str(record.seq)
            if seq not in seqs:
                seqs.add(seq)
                uniqued[rid] = record
        return uniqued

    def bind2file(self, path: FilePathType, mode: str = 'w', **kwargs):
        """
        Bind the fasta object to a file.

        But the deleted sequence will not be saved to the file.
        """
        self._binded_file = ensure_path(path).open(mode)
        self._binded_kwargs = kwargs
        self.to_fasta(self._binded_file, **kwargs)

    def unbind(self):
        """
        Unbind the fasta object from the file.
        """
        if hasattr(self, '_binded_file'):
            self._binded_file.close()
            del self._binded_file
            del self._binded_kwargs

    def __del__(self):
        self.unbind()


def read_fasta(path: FilePathType) -> Fasta:
    """
    Read fasta file as `Fasta` object.
    
    Parameters
    ----------
    path : str or Path
        The path of the fasta file.

    Returns
    ----------
    Fasta object.
    """
    if not isinstance(path, Path):
        path = Path(path)
    path = path.resolve().expanduser().absolute()
    return Fasta(map(lambda x: (x.id, x), SeqIO.parse(path, 'fasta')))


def read_seqres(path: Path) -> Fasta:
    if not isinstance(path, Path):
        path = Path(path)
    if path.suffix.lower() == '.pdb':
        seqres = SeqIO.parse(path, 'pdb-seqres')
    elif path.suffix.lower() == '.cif':
        seqres = SeqIO.parse(path, 'cif-seqres')
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return Fasta((record.id[-1], str(record.seq)) for record in seqres)


def save_fasta(
        sequences: Iterable[SeqRecord],
        path: FilePathOrIOType,
        mkdir: bool = False,
        two_line_mode: bool = False):
    """
    Save fasta file.

    Parameters
    ----------
    sequences : Iterable[SeqRecord]
        The sequences to be saved.
    path : str or Path
        The path of the fasta file to be saved.
    mkdir : bool, optional
        Whether to make the directory. The default is False.
    """
    if mkdir and isinstance(path, Path):
        path.parent.mkdir(parents=True, exist_ok=True)
    f, need_close = ensure_fileio(path, 'w')
    SeqIO.write(sequences, f, 'fasta-2line' if two_line_mode else 'fasta')
    if need_close:
        f.close()
    else:
        f.flush()


def df2fasta(df:pd.DataFrame,
             fasta_path: FilePathType, 
             *seq_cols: str, 
             id_col: Optional[str] = None, 
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


def temp_fasta(path: FilePathType, id_prefix: str = ''):
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
    subparsers = parser.add_subparsers(dest="cmd")
    csv2fasta_parser = subparsers.add_parser(
        'csv2fasta',
        help='Convert sequence csv to fasta format.')
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
    
    fasta2csv_parser = subparsers.add_parser(
        'fasta2csv',
        help='Convert sequence fasta format to csv format.')
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

    complex_parser = subparsers.add_parser(
        'complex',
        help='Create complex sequences from two fasta files.'
        )
    complex_parser.add_argument(
        '--seqs1',
        '-i1', type=Path,
        required=True,
        help='Input fasta 1')
    complex_parser.add_argument(
        '--seqs2',
        '-i2',
        type=Path,
        required=True,
        help='Input fasta 2')
    complex_parser.add_argument(
        '--output',
        '-o',
        type=Path,
        required=True,
        help='Output fasta')
    complex_parser.add_argument(
        '--linker',
        '-l',
        default='',
        help='Linker between two sequences in a complex')

    unique_parser = subparsers.add_parser(
        'unique',
        help='Remove duplicated sequences')
    unique_parser.add_argument(
        '--input',
        '-i',
        required=True,
        type=Path,
        help='Input fasta')
    unique_parser.add_argument(
        '--output',
        '-o',
        required=True,
        type=Path,
        help='Output fasta')

    args = parser.parse_args()

    if args.cmd == 'csv2fasta':
        df = pd.read_csv(args.input)
        df2fasta(df, args.output, *args.seq_cols, 
                 id_col=args.id_col,
                 mode=args.mode, sep=args.sep)

    elif args.cmd == 'fasta2csv':
        fasta = read_fasta(args.input)
        fasta.to_csv(args.output)

    elif args.cmd == 'complex':
        seqs1 = read_fasta(args.seqs1)
        seqs2 = read_fasta(args.seqs2)

        save_fasta(
            cross_create(seqs1.values(), seqs2.values(), args.linker),
            args.output)

    elif args.cmd == 'unique':
        seqs = read_fasta(args.input)
        seqs_uniqued = seqs.unique()
        print(f'Input size: {len(seqs)}, output size: {len(seqs_uniqued)}')
        seqs_uniqued.to_fasta(args.output)

    elif args.cmd is None:
        parser.print_help()

    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")
