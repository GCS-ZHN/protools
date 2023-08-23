import pandas as pd

from typing import overload, Dict, Iterable, Tuple, Union
from pathlib import Path

FilePath = Union[str, Path]


@overload
def save_fasta(seq_dict: Dict[str, str], fasta_path: FilePath):
    ...


@overload
def save_fasta(seq_iter: Iterable[Tuple[str, str]], fasta_path: FilePath):
    ...


def iter_fasta(fasta_path: FilePath) -> Iterable[Dict]:
    ...


def read_fasta(fasta_path: FilePath, id_as_index: bool = False) -> pd.DataFrame:
    ...


def df2fasta(df:pd.DataFrame,
             fasta_path: FilePath,
             *seq_cols: str,  
             id_col: str = None,
             mode: str = 'seperate', 
             sep: str = ''):
    ...
