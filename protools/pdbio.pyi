import pandas as pd

from pathlib import Path
from typing import Callable, overload, Union, Tuple, Iterable
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity


@overload
def save_to_pdb(output_path: str, *entities: Model, remarks: Iterable[str] = None) -> None:
    ...


@overload
def save_to_pdb(output_path: str, *entities: Chain, remarks: Iterable[str] = None) -> None:
    ...


@overload
def save_to_pdb(output_path: str, *entities: Residue, remarks: Iterable[str] = None) -> None:
    ...


def download_PDB(pdb_id: str, target_path: str, server: str = 'http://ftp.wwpdb.org') -> Path:
    ...


async def async_download_PDB(pdb_id: str, target_path: str, callback: Callable) -> Path:
    ...


def read_pdb_seq(pdb_file: str) -> Iterable[Tuple[str, str, str]]:
    ...


def pdb2fasta(
        fasta_file: str,
        *pdb_files: str,
        multimer_mode: str = 'joint',
        joint_sep: str = ':') -> None:
    ...


def pdb2df(entity: Union[Structure, Model, Chain, Residue], *extra_attrs: str) -> pd.DataFrame:
    ...


def read_residue(pdb: Union[Path, str, Entity], mode='centroid') -> pd.DataFrame:
    ...
