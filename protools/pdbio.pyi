from typing import overload
from Bio.PDB.Residue import Residue, Iterable
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model


@overload
def save_to_pdb(output_path: str, *entities: Model, remarks: Iterable[str] = None) -> None:
    ...


@overload
def save_to_pdb(output_path: str, *entities: Chain, remarks: Iterable[str] = None) -> None:
    ...


@overload
def save_to_pdb(output_path: str, *entities: Residue, remarks: Iterable[str] = None) -> None:
    ...