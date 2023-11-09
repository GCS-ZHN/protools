import csv
import logging
from multiprocessing import Pool
from pathlib import Path
from itertools import product, starmap
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBParser
from .utils import ensure_path, FilePath


loger = logging.getLogger(__name__)


def calc_sasa(pdb_file: FilePath, output_dir: FilePath, model_idx: int = 0):
    loger.info(f"processing {pdb_file}")
    pdb_file = ensure_path(pdb_file)
    output_dir = ensure_path(output_dir)
    loger.debug(f"creating output_dir: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{pdb_file.stem}.csv"
    loger.debug(f"reading pdb file")
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("pdb", pdb_file)
    loger.debug(f"calculating sasa")
    sasa_calculator = ShrakeRupley()
    model = struct[model_idx]
    sasa_calculator.compute(model, level="R")
    
    loger.info(f"writing result to {output_file}")
    with open(output_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["chain", "residue", "sasa"])
        for chain in model:
            for res in chain:
                writer.writerow([chain.id, res.get_id()[1], res.sasa])


def calc_sasa_from_list(
        pdb_dir: FilePath, 
        output_dir: FilePath, 
        model_idx: int = 0, 
        num_worker: int = 1):
    pdb_dir = ensure_path(pdb_dir)
    output_dir = ensure_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    schedules = product(pdb_dir.glob('*.pdb'), [output_dir], [model_idx])

    if num_worker == 1:
        list(starmap(calc_sasa, schedules))

    elif num_worker > 1:
        with Pool(num_worker) as p:
            p.starmap(calc_sasa, schedules)

    else:
        raise ValueError("num_worker must be greater than 0")
    

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--pdb_dir", "-i", type=Path, required=True)
    parser.add_argument("--output_dir", "-o", type=Path, required=True)
    parser.add_argument("--model_idx", type=int, default=0)
    parser.add_argument("--num_worker", "-n", type=int, default=1)
    parser.add_argument("--log_level", "-l", type=str, default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=logging.getLevelName(args.log_level),
                        format='%(process)d-%(asctime)s-%(levelname)s-%(message)s')
    
    calc_sasa_from_list(args.pdb_dir, args.output_dir, args.model_idx, args.num_worker)
    loger.info("done")
