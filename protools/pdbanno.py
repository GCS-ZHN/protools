import pandas as pd

from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData
from pathlib import Path
from scipy.spatial import distance
from typing import Union

from . import pdbio


def neighbor_water_count(pdbfile: Union[Path, str], level: str = 'residue', threshold: float = 2.8) -> pd.Series:
    """Count the number of water molecules in a PDB file.

    Parameters
    ----------
    pdbfile : Union[Path, str]
        Path to the PDB file.
    level : str, optional
        Level of the annotation, by default 'residue'.
    threshold : float, optional
        Distance threshold to define a water molecule,
        by default 2.8.
    
    Returns
    ----------
    pd.Series
        Number of water molecules for each model/chain/residue.

    Examples
    ----------
    >>> from protools import pdbanno
    >>> pdbanno.neighbor_water_count('3inj.pdb')
    index
    (0, A, 10, ALA)     0
    (0, A, 100, THR)    1
    (0, A, 101, TYR)    1
    (0, A, 102, LEU)    0
    (0, A, 103, ALA)    0
                    ..
    (0, H, 95, ILE)     0
    (0, H, 96, GLU)     0
    (0, H, 97, ARG)     0
    (0, H, 98, ASP)     0
    (0, H, 99, ARG)     0
    Length: 3953, dtype: int64
    """
    parser = PDBParser(QUIET=True)
    pdbfile = Path(pdbfile).absolute().resolve().expanduser()
    structure = parser.get_structure("pdb", pdbfile)
    df = pdbio.pdb2df(structure)

    df['is_standard_res'] = df['resn'].apply(
        lambda x: x.capitalize() in IUPACData.protein_letters_3to1)
    df['is_heavy_atom'] = df['element'].apply(
        lambda x: x in ['N', 'O', 'S', 'P'])
    df['is_water'] = df['resn'].apply(lambda x: x == 'HOH')

    res_hatom = df[df['is_standard_res'] & df['is_heavy_atom']]
    res_hatom_coord = res_hatom[['x', 'y', 'z']]
    water_oxygen_coord = df[df['is_water'] & df['is_heavy_atom']][['x', 'y', 'z']]

    if level == 'model':
        index = res_hatom['model']
    elif level == 'chain':
        index = res_hatom[['model', 'chain']]
    elif level == 'residue':
        index = res_hatom[['model', 'chain', 'resi', 'resn']]
    else:
        raise ValueError(f"Unknown level: {level}")

    dist = distance.cdist(res_hatom_coord, water_oxygen_coord, 'euclidean') <= threshold
    dist = pd.DataFrame(index=index, data=dist)
    dist.index.name = 'index'
    return dist.groupby('index').any().sum(axis=1)
