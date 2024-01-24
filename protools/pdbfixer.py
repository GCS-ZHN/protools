import warnings
from typing import Optional

from Bio import BiopythonWarning

from .pdbio import save_pdb, get_structure
from .utils import ensure_path

warnings.simplefilter('ignore', BiopythonWarning)


def renumber_residue(
        pdb_file: str, 
        out_file: Optional[str] = None,
        model_idx: int = 0, 
        chain_order: Optional[list] = None, 
        start: int = 1):
    """
    Renumber the residue index (resi) of a pdb file
    with continuous numbers (start from `start`)
    according to the given order of the chains.

    Parameters
    ----------
    pdb_file : str
        The path of the pdb file to be renumbered.
    out_file : str, optional
        The path of the output file.
        If not given, the output file will be
        named as `pdb_file` with suffix `_renumbered`.
    model_idx : int, optional
        The index of the model to be renumbered.
        The default is 0. Output file will only
        contain this model.
    chain_order : list, optional
        The order of the chains.
        If not given, the order of the chains
        in the original file will be used.
        If less chains are given, the rest
        chains will be appended to the end.
    start : int, optional
        The start number of the renumbering.
        The default is 1.

    Raises
    ----------
    FileNotFoundError
        If the pdb file does not exist.

    Notes
    ----------
    Renumbering will remove the insertion code of the residues!
    """
    pdb_file = ensure_path(pdb_file)
    structure = get_structure(pdb_file)
    model = structure[model_idx]

    old_chain_order = [chain.id for chain in model]

    if chain_order is None:
        chain_order = old_chain_order
    else:
        # select chain will insert first
        for chain_id in old_chain_order:
            if chain_id not in chain_order:
                chain_order.append(chain_id)
    assert len(chain_order) == len(old_chain_order), "The chain order should contain all the chains in the model"
    
    i = 0
    for chain_id in chain_order:
        for residue in model[chain_id]:
            old_id = residue.id
            residue.id = (' ', start + i, 'TMP')
            i += 1
    
    # loop twice to avoid duplicate residue id
    for chain_id in chain_order:
        for residue in model[chain_id]:
            old_id = residue.id
            residue.id = (' ', old_id[1], ' ')

    if out_file is None:
        out_file = pdb_file.parent / f"{pdb_file.stem}_renumbered{pdb_file.suffix}" 
    
    save_pdb(out_file, *[model[chain_id] for chain_id in chain_order], remarks=[f"Renumbered from {pdb_file}"])


if __name__ == '__main__':

    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd')

    # renumber residue
    parser_renres = subparsers.add_parser('renres')
    parser_renres.add_argument('pdb_file', type=str, help='the pdb file to be renumbered')
    parser_renres.add_argument('-o', '--out_file', type=str, help='the output file')
    parser_renres.add_argument('-m', '--model', type=int, default=0, help='the model to be renumbered')
    parser_renres.add_argument('-c', '--chain_order', type=str, nargs='+', help='the order of the chains')
    parser_renres.add_argument('-s', '--start', type=int, default=1, help='the start number of the renumbering')

    args = parser.parse_args()

    if args.cmd == 'renres':
        renumber_residue(args.pdb_file, args.out_file, args.model, args.chain_order, args.start)
    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")
