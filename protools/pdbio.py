#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging

from typing import Callable, Iterable, Union, overload
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBList
from pathlib import Path

@overload
def save_to_pdb(output_path: str, *entities: Model, remarks: Iterable[str]=None) -> None:
    ...

@overload
def save_to_pdb(output_path: str, *entities: Chain, remarks: Iterable[str]=None) -> None:
    ...

@overload
def save_to_pdb(output_path: str, *entities: Residue, remarks: Iterable[str]=None) -> None:
    ...

def save_to_pdb(output_path: str, *entities: Union[Model, Chain, Residue], remarks: Iterable[str]=None) -> None:
    """
    Save entities to a PDB file.

    Parameters
    ----------
    output_path : str
        Path to the output PDB file.
    entities : Model, Chain, or Residue
        Entities to save.
    remarks : Iterable[str], optional
        Remarks to be written to the PDB file.

    Raises
    ------
    ValueError
        If no entities are provided or if the
        entities are not of the same type.

    TypeError
        If the entities are not of type Model,
        Chain, or Residue.
    """
    if len(entities) == 0:
        raise ValueError("No entities to save")
    
    def _loop_type_check(entities, *allowed_types, allow_none=False):
        for entity in entities:
            if not isinstance(entity, allowed_types):
                if not allow_none or entity is not None:
                    raise ValueError(f"Unsupported type {type(entity)}")
            yield entity
    
    if isinstance(entities[0], Model):
        structure = Structure("pdb")
        for model in _loop_type_check(entities, Model):
            structure.add(model)
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        with open(output_path, "w") as fp:
            fp.write("REMARK 220 Generated by Python\n")
            if remarks is not None:
                for remark in remarks:
                    fp.write(f"REMARK 999 {remark}\n")
            pdb_io.save(fp)

    elif isinstance(entities[0], Chain):
        models = [Model("model_0")]
        for chain in _loop_type_check(entities, Chain, allow_none=True):
            if chain is None:
                models.append(Model(f"model_{len(models)}"))
            else:
                models[-1].add(chain)
        save_to_pdb(output_path, *models, remarks=remarks)

    elif isinstance(entities[0], Residue):
        chains = [Chain("A")]
        for residue in _loop_type_check(entities, Residue, allow_none=True):
            if residue is None:
                chain_id = chr(ord(chains[-1].get_id()) + 1)
                if chain_id > "Z":
                    raise ValueError("Too many chains")
                chains.append(Chain(chain_id))
            else:
                chains[-1].add(residue)
        save_to_pdb(output_path, *chains, remarks=remarks)

    else:
        raise TypeError(f"Unsupported type {type(entities[0])}")


def download_PDB(pdb_id: str, target_path: str, server: str='http://ftp.wwpdb.org') -> Path:
    """
    Download a PDB file from the PDB database.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the PDB file to download.

    target_path : str
        Path to save the downloaded PDB file.
        If current pdb file exists, download will be skipped.

    server : str, optional
        URL of the PDB server to download from.
        Default is 'http://ftp.wwpdb.org'.

    Returns
    ----------
    file : Path
        Path to the downloaded PDB file.

    Raises
    ----------
    FileNotFoundError
        If the PDB file could not be downloaded.
    """
    pdbl = PDBList(server=server, verbose=False)
    file = pdbl.retrieve_pdb_file(pdb_id, pdir=target_path, file_format='pdb')
    file = Path(file)
    if not file.is_file():
        raise FileNotFoundError(f"Could not download PDB {pdb_id} to {target_path}")
    return file


async def async_download_PDB(pdb_id: str, target_path: str, callback: Callable) -> Path:
    """
    Asynchronously download a PDB file from the PDB database.
    Downloading is a IO-bound task, so it is suitable to be
    run asynchronously.
    
    Parameters
    ----------
    pdb_id : str
        PDB ID of the PDB file to download.

    target_path : str
        Path to save the downloaded PDB file.
        If current pdb file exists, download will be skipped.

    callback : Callable
        Callback function to be called after the download
        is finished, no arguments are passed to the callback;
        if the download failed, the PDB ID will be passed.

    Returns
    ----------
    file : Path
        Path to the downloaded PDB file.

    Notes
    ----------
    This function is not thread-safe, so it should be run
    in a single thread. Downloading failed will not raise
    an exception, which is different from the method 
    `download_PDB`.
    """
    try:
        file = download_PDB(pdb_id, target_path)
        callback()
        return file
    except RuntimeError as e:
        callback(pdb_id)
        logging.error(e)


if __name__ == "__main__":
    import asyncio
    from tqdm.auto import tqdm
    from argparse import ArgumentParser
    parser = ArgumentParser()

    # add subcommand 'download'
    subparsers = parser.add_subparsers(dest="subcommand")
    download_parser = subparsers.add_parser("download")

    pdb_id_group = download_parser.add_mutually_exclusive_group(required=True)
    pdb_id_group.add_argument(
        "--pdb_ids",
        "-i",
        nargs="+",
        help="PDB IDs to download")
    pdb_id_group.add_argument(
        "--pdb_id_file",
        "-f",
        type=str,
        help="Path to a file containing PDB IDs to download")

    download_parser.add_argument(
        "--target_path",
        "-d",
        type=str,
        default=".",
        help="Path to save the downloaded PDB files")

    # add subcommand 'help'
    help_parser = subparsers.add_parser("help")
    help_parser.add_argument("command", type=str, help="Subcommand to show help for")

    args = parser.parse_args()

    if args.subcommand == "download":
        if args.pdb_ids is not None:
            pdb_ids = args.pdb_ids
        else:
            with open(args.pdb_id_file, "r") as fp:
                pdb_ids = [line.strip() for line in fp]
        
        for pid in pdb_ids:
            if len(pid)!=4:
                raise ValueError(f"Invalid pdbid {pid}")

        output_path = Path(args.target_path).resolve()
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        elif not output_path.is_dir():
            raise FileNotFoundError(f"Target path {args.target_path} is not a directory")

        process_bar = tqdm(total=len(pdb_ids), desc="Downloading PDB files")
        download_failed = []
        def callback(*args):
            if len(args) > 0:
                download_failed.append(args[0])
            process_bar.update(1)

        tasks = [async_download_PDB(pdb_id, str(output_path), callback) for pdb_id in pdb_ids]
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.gather(*tasks))
        process_bar.close()
    elif args.subcommand == "help":
        subparsers.choices[args.command].print_help()
    else:
        parser.print_help()
