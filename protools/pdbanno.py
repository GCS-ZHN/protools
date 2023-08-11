import pandas as pd

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import IUPACData
from pathlib import Path
from scipy.spatial import distance
from typing import Union

from . import pdbio


def neighbor_water_count(
        pdb: Union[Path, str, Structure], 
        level: str = 'R', 
        threshold: float = 2.8) -> Union[int, pd.Series]:
    """Count the number of water molecules in a PDB file.

    Parameters
    ----------
    pdbfile : Union[Path, str, Structure]
        Path to the PDB file or a Bio.PDB.Structure.Structure object.
    level : str, optional
        Level of the annotation, by default 'R'.
        'S' for structure, 'M' for model, 'C' for chain, 
        'R' for residue.
    threshold : float, optional
        Distance threshold to define a water molecule,
        by default 2.8.
    
    Returns
    ----------
    pd.Series or int
        Number of water molecules for each model/chain/residue.
        If `level` is 'S', return the total number of water molecules.

    Examples
    ----------
    >>> from protools import pdbanno
    >>> pdbanno.neighbor_water_count('3inj.pdb')
    model  chain  resi  resn
    0      A      10    ALA     0
                  100   THR     1
                  101   TYR     1
                  102   LEU     0
                  103   ALA     0
                              ..
           H      95    ILE     0
                  96    GLU     0
                  97    ARG     0
                  98    ASP     0
                  99    ARG     0
    Length: 3953, dtype: int64
    """
    parser = PDBParser(QUIET=True)
    if isinstance(pdb, Structure):
        structure = pdb
    elif isinstance(pdb, (str, Path)):
        pdb = Path(pdb).absolute().resolve().expanduser()
        structure = parser.get_structure("pdb", pdb)
    else:
        raise TypeError(f"Unknown type: {type(pdb)}")

    df = pdbio.pdb2df(structure)

    df['is_standard_res'] = df['resn'].apply(
        lambda x: x.capitalize() in IUPACData.protein_letters_3to1)
    df['is_heavy_atom'] = df['element'].apply(
        lambda x: x in ['N', 'O', 'S', 'P'])
    df['is_water'] = df['resn'].apply(lambda x: x == 'HOH')

    res_hatom = df[df['is_standard_res'] & df['is_heavy_atom']]
    res_hatom_coord = res_hatom[['x', 'y', 'z']]
    water_oxygen_coord = df[df['is_water'] & df['is_heavy_atom']][['x', 'y', 'z']]

    if level in ('S', 'M'):
        level_columns = ['model']
    elif level == 'C':
        level_columns = ['model', 'chain']
    elif level == 'R':
        level_columns = ['model', 'chain', 'resi', 'resn']
    else:
        raise ValueError(f"Unknown level: {level}")


    dist = distance.cdist(res_hatom_coord, water_oxygen_coord, 'euclidean') <= threshold
    dist = pd.DataFrame(index=res_hatom_coord.index, data=dist)
    if level == 'S':
        return dist.any().sum()

    dist = dist.join(res_hatom[level_columns])
    return dist.groupby(level_columns).any().sum(axis=1)


def calc_sasa(pdb: Structure, radius: float = 1.4, standard: bool = True) -> pd.Series:
    """
    calculate the solvent accessible surface area (SASA) of each residue
    in a PDB file.

    Parameters
    ----------
    pdb : Structure
        The PDB structure.
    radius : float, optional
        The probe radius, by default 1.4.
    standard : bool, optional
        Only calculate the SASA of standard amino acids,
        by default True.

    Returns
    ----------
    pd.Series
        The SASA of each residue.
        The index is a tuple of (model, chain, resi, resn).
    """
    sasa_calc = ShrakeRupley(probe_radius=radius)
    sasa_calc.compute(pdb, level='R')
    result = dict()
    for model in pdb:
        for chain in model:
            for res in chain:
                if standard and res.resname.capitalize() not in IUPACData.protein_letters_3to1:
                    continue

                if res.resname == 'HOH':
                    continue

                resid = f'{res.id[1]}{res.id[2]}'.strip()
                index = (model.id, chain.id, resid, res.resname)
                result[index] = res.sasa
    return pd.Series(result)
