#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Tuple
from scipy.spatial.distance import cdist

from protools.aa import RESIDUE2SIDECHAIN_3LETTER
from protools.pdbio import Structure, get_structure, save_pdb, pdb2df
from protools.typedef import FilePathType, StructureFragmentAAType
from protools.utils import ensure_path


def parser_resi(resi: str) -> Iterable[Tuple]:
    """
    Parse a string of residue selection into a list of residue range.

    Parameters
    ----------
    resi : str
        A string of residue selection. different range
        should be separated by comma. Each range should
        be start with a chain id.


    Examples
    ----------
    >>> list(parser_resi('H1A-4,F6,M8-10'))
    [('H', ('1A', '4')), ('F', ('6', '6')), ('M', ('8', '10'))]
    """
    for resi_range in resi.split(','):
        chain_id = resi_range[0]
        resi_range = resi_range[1:]
        start, *end = resi_range.split('-')
        if len(end) == 0:
            end = start
        elif len(end) == 1:
            end = end[0]
        else:
            raise ValueError(f"Invalid residue range: {resi_range}")
        yield chain_id, (start, end)


def _build_residue_dict(
        structure: Structure, 
        model_id: int = 0) -> Tuple[dict, dict]:
    """
    Internal function to build a dict of residue index to residue object
    and a dict of residue object to residue index.
    """
    resi2idx = {}
    idx2res = []
    model = structure[model_id]
    for idx, residue in enumerate(model.get_residues()):
        chain_id = residue.get_parent().get_id()
        resi = (chain_id, ''.join(map(str, residue.get_id())).strip())
        resi2idx[resi] = idx
        idx2res.append(residue)
    return resi2idx, idx2res


def extract(
        pdb_file: FilePathType, 
        out_file: FilePathType,
        resi_selection: str, 
        remain: bool = True,
        model_id: int = 0):
    """
    Extract residues from pdb file.

    Parameters
    ----------
    pdb_file : FilePathType
        Path to the pdb file.
    out_file : FilePathType
        Path to the output pdb file.
    resi_selection : str
        A string of residue selection. different range
        should be separated by comma. Each range should
        be start with a chain id.
    remain : bool, optional
        Whether to remain the selected residues or not,
        by default True.
    model_id : int, optional
        The index of the model to be used, by default 0.
    """
    pdb_file = ensure_path(pdb_file)
    out_file = ensure_path(out_file)
    structure = get_structure(pdb_file)
    resi2idx, idx2res = _build_residue_dict(structure, model_id)
    idx_mask = [not remain] * len(idx2res)
    for chain_id, (start, end) in parser_resi(resi_selection):
        start = resi2idx[(chain_id, start)]
        end = resi2idx[(chain_id, end)]
        for idx in range(start, end + 1):
            idx_mask[idx] = remain

    chains = []
    for idx, residue in enumerate(idx2res):
        chain = residue.get_parent()
        if idx_mask[idx]:
            if chain not in chains:
                chains.append(chain)
            continue
        chain.detach_child(residue.get_id())

    save_pdb(out_file, *chains)


def batch_extract(
        pdb_files: Iterable[FilePathType], 
        out_dir: FilePathType,
        resi_selection: str, 
        remain: bool = True,
        model_id: int = 0,
        num_process: int = 1):
    """
    Extract residues from pdb files.
    """
    pdb_files = map(ensure_path, pdb_files)
    out_dir = ensure_path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pool = Pool(num_process)

    for pdb_file in pdb_files:
        out_file = out_dir / (pdb_file.stem + '.pdb')
        pool.apply_async(
            extract,
            args=(pdb_file, out_file, resi_selection, remain, model_id)
            )
    pool.close()
    pool.join()


def _get_sidechain_radius(data: pd.DataFrame) -> pd.Series:
    """
    Get the sidechain radius of each residue in the dataframe.
    The sidechain radius is defined as the maximum distance between
    the CA atom and any sidechain atom.
    If the residue has no sidechain atoms, the radius is 0.
    If the residue has no CA atom, the radius is infinity."""
    def _get_radius(grouped_data: pd.DataFrame) -> float:
        resn = grouped_data.name[3].capitalize()
        if resn not in RESIDUE2SIDECHAIN_3LETTER:
            raise ValueError(f"Unknown residue name: {resn}")
        sidechain_atoms = RESIDUE2SIDECHAIN_3LETTER[resn]
        if len(sidechain_atoms) == 0:
            return 0.0
        sidechain_data = grouped_data[grouped_data['name'].isin(sidechain_atoms)]
        ca_data = grouped_data[grouped_data['name'] == 'CA']
        if len(ca_data) == 0 or len(sidechain_data) == 0:
            return float('inf')
        ca_coord = ca_data[['x', 'y', 'z']].to_numpy()
        sidechain_coord = sidechain_data[['x', 'y', 'z']].to_numpy()
        dists = cdist(ca_coord, sidechain_coord)
        return dists.max()
    return data.groupby(data.index).apply(_get_radius)


def distance(
        entity1: StructureFragmentAAType,
        entity2: StructureFragmentAAType,
        dist_type: str = 'ca') -> pd.DataFrame:
    entity1_df = pdb2df(entity1)
    entity2_df = pdb2df(entity2)
    entity1_df.set_index(['model', 'chain', 'seqid', 'resn'], inplace=True)
    entity2_df.set_index(['model', 'chain', 'seqid', 'resn'], inplace=True)

    if dist_type == 'ca':
        entity1_ca_df = entity1_df[entity1_df['name'] == 'CA']
        entity2_ca_df = entity2_df[entity2_df['name'] == 'CA']
        dist = cdist(
            entity1_ca_df[['x', 'y', 'z']].to_numpy(),
            entity2_ca_df[['x', 'y', 'z']].to_numpy())
        dist_df = pd.DataFrame(
            dist,
            index=entity1_ca_df.index,
            columns=entity2_ca_df.index
            )
    elif dist_type == 'full_atom':
        dist = cdist(
            entity1_df[['x', 'y', 'z']].to_numpy(),
            entity2_df[['x', 'y', 'z']].to_numpy())
        dist_df = pd.DataFrame(
            dist,
            index=entity1_df.index,
            columns=entity2_df.index)
        dist_df = dist_df.groupby(dist_df.index).min().T.groupby(dist_df.columns).min().T
        dist_df.index = pd.MultiIndex.from_tuples(
            dist_df.index, names=entity1_df.index.names)
        dist_df.columns = pd.MultiIndex.from_tuples(
            dist_df.columns, names=entity2_df.index.names)
    elif dist_type == 'sidechain_radius':
        entity1_ca_df = entity1_df[entity1_df['name'] == 'CA']
        entity2_ca_df = entity2_df[entity2_df['name'] == 'CA']
        entity1_sidechain_radius = _get_sidechain_radius(entity1_df)
        entity2_sidechain_radius = _get_sidechain_radius(entity2_df)
        dist = cdist(
            entity1_ca_df[['x', 'y', 'z']].to_numpy(),
            entity2_ca_df[['x', 'y', 'z']].to_numpy())
        dist -= entity1_sidechain_radius.values[:, None]
        dist -= entity2_sidechain_radius.values[None, :]
        dist[dist < 0] = 0
        dist_df = pd.DataFrame(
            dist,
            index=entity1_ca_df.index,
            columns=entity2_ca_df.index)
    else:
        raise ValueError(f"Unknown distance type: {dist_type}")
    return dist_df


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(
        '--pdb-file',
        '-i',
        type=Path,
        nargs='+',
        required=True
        )
    
    parser.add_argument(
        '--out-dir',
        '-o',
        type=Path,
        required=True
        )
    
    parser.add_argument(
        '--resi-selection',
        '-r',
        type=str,
        required=True
        )
    
    parser.add_argument(
        '--remain',
        '-m',
        action='store_true'
        )
    
    parser.add_argument(
        '--model-id',
        '-d',
        type=int,
        default=0
        )
    
    parser.add_argument(
        '--num-process',
        '-p',
        type=int,
        default=1
        )
    
    args = parser.parse_args()
    batch_extract(
        args.pdb_file,
        args.out_dir,
        args.resi_selection,
        args.remain,
        args.model_id,
        args.num_process
        )
