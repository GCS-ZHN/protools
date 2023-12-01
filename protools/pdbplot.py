from .utils import require_package
from pathlib import Path
from typing import Optional
import logging

try:
    require_pymol = require_package("pymol", "conda install -c schrodinger pymol")
    from pymol import cmd
except ImportError:
    pass

logger = logging.getLogger(__name__)


@require_pymol
def align_all_designs(
    pdb_dir: Path, 
    output_name: Path,
    receptor_chain: Optional[str] = None, 
    ref_pdb: Optional[Path] = None):

    pdbs = pdb_dir.glob("*.pdb")
    remove_ref = True
    if ref_pdb is None:
        ref_pdb = next(pdbs)
        remove_ref = False
    if not ref_pdb.is_file():
        raise FileNotFoundError(f"{ref_pdb} does not exist.")
    
    cmd.load(str(ref_pdb), ref_pdb.stem)
    logger.info(f'loaded {ref_pdb.stem} as reference.')

    for pdb in pdbs:
        cmd.load(str(pdb), pdb.stem)
        logger.info(f'loaded {pdb.stem} as design.')
        align_from_str = pdb.stem
        align_to_str = ref_pdb.stem

        if receptor_chain is not None:
            align_from_str += f" and chain {receptor_chain}"
            align_to_str += f" and chain {receptor_chain}"

        cmd.align(align_from_str, align_to_str)
    
    if remove_ref:
        cmd.delete(ref_pdb.stem)

    output_name.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f'saving aligned structures to {output_name}')
    cmd.save(str(output_name))
    cmd.delete("all")


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")
    design_align_parser = subparsers.add_parser("align")
    design_align_parser.add_argument("--pdb_dir" ,'-i', type=Path, required=True)
    design_align_parser.add_argument("--output_name" , '-o', type=Path, required=True)
    design_align_parser.add_argument("--receptor_chain", '-r', type=str)
    design_align_parser.add_argument("--ref_pdb", '-p', type=Path)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, 
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    if args.cmd == "align":
        align_all_designs(
            args.pdb_dir, 
            args.output_name, 
            receptor_chain=args.receptor_chain,
            ref_pdb=args.ref_pdb
        )

    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")
