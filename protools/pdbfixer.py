from Bio.PDB import PDBParser
from pathlib import Path
from .utils import FilePath, ensure_path
from .pdbio import save_to_pdb


def reset_resi(pdb_path: FilePath):
    parser = PDBParser(QUIET=True)
    pdb_path = ensure_path(pdb_path)
    structure = parser.get_structure("pdb", pdb_path)
    for model in structure:
        for i, residue in enumerate(model.get_residues(), 1):
            residue.id = (' ', i, ' ')
    save_to_pdb(
        pdb_path.with_name(pdb_path.stem + "-fixed.pdb"),
        structure
    )


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("data_path", type=Path)
    args = parser.parse_args()

    reset_resi(args.data_path)
