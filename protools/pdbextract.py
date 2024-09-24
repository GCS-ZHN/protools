#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Tuple

from .pdbio import Structure, get_structure, save_pdb
from .typedef import FilePathType
from .utils import ensure_path


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
