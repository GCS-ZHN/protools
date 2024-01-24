from pathlib import Path
from typing import Union
from io import IOBase
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Type definitions
FilePathType = Union[Path, str]
FilePathOrIOType = Union[FilePathType, IOBase]
StructureFragmentType = Union[Structure, Model, Chain]
StructureFragmentAAType =Union[Structure, Model, Chain, Residue]
SeqLikeType = Union[str, Seq, SeqRecord]
