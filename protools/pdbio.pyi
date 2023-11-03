import pandas as pd

from pathlib import Path
from typing import Callable, overload, Union, Tuple, Iterable
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
from .utils import FilePath


@overload
def save_to_pdb(output_path: FilePath, *entities: Structure, remarks: Iterable[str] = None) -> None:
    ...


def save_to_pdb(output_path: FilePath, *entities: Model, remarks: Iterable[str] = None) -> None:
    ...


@overload
def save_to_pdb(output_path: FilePath, *entities: Chain, remarks: Iterable[str] = None) -> None:
    ...


@overload
def save_to_pdb(output_path: FilePath, *entities: Residue, remarks: Iterable[str] = None) -> None:
    ...


def download_PDB(pdb_id: str, target_path: FilePath, server: str = 'http://ftp.wwpdb.org') -> Path:
    ...


async def async_download_PDB(pdb_id: str, target_path: FilePath, callback: Callable) -> Path:
    ...


def read_pdb_seq(pdb_file: FilePath) -> Iterable[Tuple[str, str, str]]:
    ...


def pdb2fasta(
        fasta_file: FilePath,
        *pdb_files: FilePath,
        multimer_mode: str = 'joint',
        joint_sep: str = ':') -> None:
    ...


def pdb2df(entity: Union[Structure, Model, Chain, Residue], *extra_attrs: str) -> pd.DataFrame:
    ...


def read_residue(pdb: Union[Path, str, Entity], mode='centroid') -> pd.DataFrame:
    ...


def get_pdb_remarks(pdb_file: FilePath) -> Iterable[str]:
    ...