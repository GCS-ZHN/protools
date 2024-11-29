import pandas as pd
import tempfile
import csv

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
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
            if hasattr(self, '_binded_fileio'):
                save_fasta(
                    [__value],
                    self._binded_fileio, 
                    mode=self._binded_fileio.mode,
                    **self._binded_kwargs)
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
    
    def to_fasta(self, path: FilePathOrIOType,
                 mkdir: bool = False, mode: str = 'w', **kwargs):
        """
        Save as FASTA format.
        """
        save_fasta(self.values(), 
            path=path, mkdir=mkdir, mode=mode, **kwargs)

    def to_csv(self, path: FilePathOrIOType,
               mkdir: bool = False, mode: str = 'w'):
        """
        Save as CSV format.
        """
        f, need_close = ensure_fileio(path, mode)
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

    def bind2file(self, path: FilePathType, **kwargs):
        """
        Bind the fasta object to a file.

        But the deleted sequence will not be saved to the file.
        """
        _binded_fileio = ensure_path(path).open('a+')
        if _binded_fileio.tell() > 0:
            if len(self) > 0:
                raise ValueError('Both file and fasta object are not empty.')
            _binded_fileio.seek(0)
            self.update(read_fasta(_binded_fileio, 'a+'))
        if len(self) > 0:
            # refresh file old content with current kwargs
            _binded_fileio = ensure_path(path).open('w+')
            self.to_fasta(
                _binded_fileio, mode='w+', **kwargs)
        self._binded_fileio = _binded_fileio
        self._binded_kwargs = kwargs

    def unbind(self):
        """
        Unbind the fasta object from the file.
        """
        if hasattr(self, '_binded_file'):
            self._binded_fileio.close()
            del self._binded_fileio
            del self._binded_kwargs

    def __del__(self):
        self.unbind()


def read_fasta(path: FilePathOrIOType, mode: str = 'r') -> Fasta:
    """
    Read fasta file as `Fasta` object.
    
    Parameters
    ----------
    path : str, Path or file-like object
        The path of the fasta file.

    Returns
    ----------
    Fasta object.
    """
    path, need_close = ensure_fileio(path, mode)
    res = Fasta(map(lambda x: (x.id, x), SeqIO.parse(path, 'fasta')))
    if need_close:
        path.close()
    return res


def read_a3m(path: FilePathOrIOType, is_multimer: bool = False) -> Fasta:
    """
    Read multiple-sequence alignment files generated by ColabSearch.

    Parameters
    ----------
    path : str, Path or file-like object
        The path of the fasta file.

    is_multimer: bool
        Whether a multimer complex.
    Returns
    ----------
    Fasta object or a list Fasta objects.
    """
    if not is_multimer:
        return read_fasta(path)
    
    path, need_close = ensure_fileio(path)
    # seqeunce size record in first line, e.g. #121,221,33
    size = path.readline()[1:].split()[0].split(',')
    size = list(map(int, size))
    res = [Fasta() for _ in size]
    query_names = None
    # only for not paired msa
    current_query = None
    for rid, heads in enumerate(map(lambda x: x.strip()[1:].split(), path)):
        # must before any continue block to be sure paired
        seq = path.readline().strip()

        if rid == 0:
            query_names = heads

        if len(heads) == 1 and heads[0] in query_names:
            current_query = heads[0]
            continue
        
        # means follow seq are note paired
        if current_query is not None:
            head = heads[0]
            heads = ['DUMMY'] * 3
            heads[query_names.index(current_query)] = head

        aligned_length = 0
        start_idx = 0
        current_no = 0
        for idx, s in enumerate(seq):
            if s.isupper() or s == '-':
                aligned_length += 1
            elif not s.islower():
                raise ValueError(f'Invalid symol for seq {heads} at position {idx+1}')
            if aligned_length == size[current_no]:
                current_seq = seq[start_idx: idx + 1]
                if heads[current_no] != 'DUMMY':
                    assert any(s!='-' for s in current_seq)
                    head = heads[current_no]
                    uid = 1
                    # add different fragments from the same database entry
                    while head in res[current_no] and \
                    str(res[current_no][head].seq) != current_seq:
                        head = f'{heads[current_no]}--unique-{uid}'
                        uid += 1
                    res[current_no][head] = current_seq
                else:
                    assert all(s=='-' for s in current_seq)
                start_idx = idx + 1
                current_no += 1
                aligned_length = 0
        assert current_no == len(size)
        assert aligned_length == 0
        assert start_idx == len(seq)
    assert current_query == query_names[-1]

    if need_close:
        path.close()
    return res

def read_pdb_seqres(path: Path) -> Dict[str, Iterable[str]]:
    """
    Read sequence from PDB SEQRES record.
    keep the format without any conversion.
    any unstandard amino acid will be kept as is.

    Parameters
    ----------
    path : Path
        The path of the PDB file.

    Returns
    ----------
    Dict[str, Iterable[str]]
        The sequence data. The key is the chain id.
        The value is the list of 3-letter amino acid 
        code or nucleotide code.
    """
    data = {}
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith('SEQRES'):
                chain_id = line[11]
                seq = line[19:].split()
                data.setdefault(chain_id, []).extend(seq)
    return data


def read_mmcif_seqres(path: Path, auth: bool = True) -> Dict[str, Iterable[str]]:
    """
    Read sequence from mmCIF file.
    Similar to `read_pdb_seqres`.

    Parameters
    ----------
    path : Path
        The path of the mmCIF file.
    auth : bool, optional
        Whether to use the auth chain id.
        The default is True.

    Returns
    ----------
    Dict[str, Iterable[str]]
        The sequence data. The key is the chain id.
        The value is the list of 3-letter amino acid 
        code or nucleotide code.

    See Also
    ----------
    read_pdb_seqres
    """
    mmcif_dict = MMCIF2Dict(path)
    if auth:
        chain_list = mmcif_dict['_pdbx_poly_seq_scheme.pdb_strand_id']
    else:
        chain_list = mmcif_dict['_pdbx_poly_seq_scheme.asym_id']
    aa_list = mmcif_dict['_pdbx_poly_seq_scheme.mon_id']
    data = {}
    for chain_id, seq in zip(chain_list, aa_list):
        data.setdefault(chain_id, []).append(seq)
    return data


def read_seqres(path: Path, auth: bool = True) -> Fasta:
    """
    Read sequence from PDB or mmCIF file.
    Extended from BioPython. Only support
    amino acid sequence.

    Parameters
    ----------
    path : Path
        The path of the PDB or mmCIF file.
    auth : bool, optional
        Whether to use the auth chain id.
        The default is True.

    Returns
    ----------
    Fasta
        The sequence data. The key is the chain id.
        The value is the sequence record.

    Notes
    ----------
    This function only support amino acid sequence, and
    return one-letter amino acid sequence as Fasta object,
    which may result in loss of information.
    """
    if not isinstance(path, Path):
        path = Path(path)
    if path.suffix.lower() == '.pdb':
        seqres = SeqIO.parse(path, 'pdb-seqres')
        seqres_iter = (
            (record.id.split(':')[-1], str(record.seq)) for record in seqres)

    elif path.suffix.lower() == '.cif':
        seqres = SeqIO.parse(path, 'cif-seqres')
        seqres_iter = (
            (record.id.split(':')[-1], str(record.seq)) for record in seqres)
        # as bioython use '_pdbx_poly_seq_scheme.asym_id' in 
        # `Bio.SeqIO.PdbIO.CifSeqresIterator`, we need to convert
        # chain id to auth chain id by '_pdbx_poly_seq_scheme.pdb_strand_id'
        if auth:
            mmcif_dict = MMCIF2Dict(path)
            label2auth = dict(zip(
                mmcif_dict['_pdbx_poly_seq_scheme.asym_id'],
                mmcif_dict['_pdbx_poly_seq_scheme.pdb_strand_id']))
            seqres_iter = (
                (label2auth[label], seq) for label, seq in seqres_iter)

    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return Fasta(seqres_iter)


def save_fasta(
        sequences: Iterable[SeqRecord],
        path: FilePathOrIOType,
        mkdir: bool = False,
        two_line_mode: bool = False,
        mode: str = 'w'):
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
    f, need_close = ensure_fileio(path, mode)
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
