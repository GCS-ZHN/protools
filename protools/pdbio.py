#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
import pandas as pd

from typing import Callable, Iterable, Tuple, Union
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Entity import Entity
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBList, PDBParser
from Bio.SeqUtils import IUPACData
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from fasta import FASTA


__all__ = [
    'save_to_pdb',
    'download_PDB',
    'async_download_PDB',
    'read_pdb_seq',
    'pdb2fasta',
    'pdb2df']


def save_to_pdb(output_path: str, *entities: Union[Model, Chain, Residue], remarks: Iterable[str] = None) -> None:
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


def download_PDB(pdb_id: str, target_path: str, server: str = 'http://ftp.wwpdb.org') -> Path:
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
        raise FileNotFoundError(
            f"Could not download PDB {pdb_id} to {target_path}")
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


def read_pdb_seq(pdb_file: str) -> Iterable[Tuple[str, str, str]]:
    """
    Extract the sequence of a PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    ----------
    seq_iter : Iterable[Tuple[str, str, str]]
        An iterable of tuples (model_id, chain_id, seq).
    """
    pdb_file = Path(pdb_file).resolve()
    if not pdb_file.exists():
        raise FileNotFoundError(f"Could not find PDB file {pdb_file}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)
    for model in structure:
        for chain in model:
            seq = "".join(IUPACData.protein_letters_3to1[residue.get_resname(
            ).capitalize()] for residue in chain)

            yield model.id, chain.id, seq


def pdb2fasta(
        fasta_file: str,
        *pdb_files: str,
        multimer_mode: str = 'joint',
        joint_sep: str = ':') -> None:
    """
    Convert PDB files to a fasta file.

    Parameters
    ----------
    fasta_file : str
        Path to the output fasta file.

    pdb_files : str
        Paths to the PDB files to be converted.
        At least one PDB file should be provided.

    multimer_mode : str, optional
        Mode of the conversion. 'seperate' for
        single sequence per entry, 'joint' for
        joint all sequences in a PDB complex file.

    joint_sep : str, optional
        Separator of the joint sequence. The default is ':'.
        Only works when multimer_mode is 'joint'.
        Different tools may have different requirements
        to process the joint sequence.

    """
    if len(pdb_files) == 0:
        raise ValueError("No PDB files to convert")

    fasta_file = Path(fasta_file).resolve()
    if not fasta_file.parent.exists():
        fasta_file.parent.mkdir(parents=True, exist_ok=True)

    with FASTA(str(fasta_file)) as f:
        for pdb_file in pdb_files:
            pdb_id = Path(pdb_file).stem
            seq_iter = read_pdb_seq(pdb_file)
            if multimer_mode == 'seperate':
                for model_id, chain_id, seq in seq_iter:
                    seq_id = f"{pdb_id}_{model_id}_{chain_id}"
                    seq_record = SeqRecord(Seq(seq), id=seq_id, description='')
                    f.add_seq(seq_record)
            elif multimer_mode == 'joint':
                seq = joint_sep.join(seq for _, _, seq in seq_iter)
                seq_id = pdb_id
                seq_record = SeqRecord(Seq(seq), id=seq_id, description='')
                f.add_seq(seq_record)
            else:
                raise ValueError(
                    f"Unsupported multimer mode {multimer_mode}")


def pdb2df(entity: Union[Structure, Model, Chain, Residue], *extra_attrs: str) -> pd.DataFrame:
    """
    Convert an entity to a pandas DataFrame.

    Parameters
    ----------
    entity : Union[Structure, Model, Chain, Residue]
        Entity to be converted.

    extra_attrs : str
        Extra atom attributes to be added to the DataFrame.
        For example, 'bfactor', 'occupancy', 'altloc', etc.

    Returns
    ----------
    df : pd.DataFrame
        DataFrame of the entity.
    """
    def _atom_to_dict(atom:Atom):
        coord = atom.get_coord()
        residue = atom.get_parent()
        resids = residue.get_id()
        chain = residue.get_parent()
        model = chain.get_parent()
        res = {
            'id': atom.get_serial_number(), 
            'name': atom.get_name(),
            'resn': residue.get_resname(),
            'resi': f'{resids[1]}{resids[2]}'.strip(),
            'chain': chain.get_id(), 
            'x': coord[0], 
            'y': coord[1], 
            'z': coord[2],
            'element': atom.element,
            'model': model.get_id()
        }
        for attr in extra_attrs:
            res[attr] = getattr(atom, attr)
        return res
    
    return pd.DataFrame(_atom_to_dict(atom) for atom in entity.get_atoms())


def read_residue(pdb: Union[Path, str, Entity], mode='centroid') -> pd.DataFrame:
    """
    Read a PDB file and return a DataFrame of the residues.

    Parameters
    ----------
    pdb_file : Union[Path, str]
        Path to the PDB file.

    mode : str, optional
        Mode of the residue coordinates. 'centroid' for
        the centroid center of the residue, 'fuc' for the
        first atom of the residue, 'CA' for the alpha carbon.

    Returns
    ----------
    df : pd.DataFrame
        DataFrame of the residues.

    Notes
    ----------
    Atoms missing in PDB will not included. So please
    fix the PDB file before using this function.

    Examples
    ----------
    >>> df = read_residue('1a12.pdb', mode='centroid')
    >>> df.head()
                                  x          y          z
    model chain resi resn                                
    0     A     100  GLY   7.494460 -33.223431  16.242475
                101  ARG   3.844297 -34.536686  13.638067
                102  ASP   2.507354 -39.138135  16.654611
                103  THR   3.519204 -40.032841  12.495594
                104  SER   1.635237 -43.545750  14.313151
    """
    if isinstance(pdb, Entity):
        structure = pdb
    elif isinstance(pdb, (Path, str)):
        pdb = Path(pdb).resolve().absolute().expanduser()
        if not pdb.exists():
            raise FileNotFoundError(f"Could not find PDB file {pdb}")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("pdb", pdb)
    else:
        raise TypeError(f"Unsupported type {type(pdb)}")

    if mode == 'centroid':
        df = pdb2df(structure, 'mass')
        df['x'] *= df['mass']
        df['y'] *= df['mass']
        df['z'] *= df['mass']
        df = df.drop(['id', 'name', 'element'], axis=1)
        df = df.groupby(['model', 'chain', 'resi', 'resn']).sum()
        df['x'] /= df['mass']
        df['y'] /= df['mass']
        df['z'] /= df['mass']
        df = df.drop('mass', axis=1)
        return df
    
    if mode in ('fuc', 'CA'):
        df = pdb2df(structure)
        if mode == 'CA':
            df = df[df['name'] == 'CA']
        df = df.drop(['id', 'name', 'element'], axis=1)
        df = df.groupby(['model', 'chain', 'resi', 'resn']).mean()
        return df
    
    raise ValueError('mode must be one of "centroid", "fuc", "CA"')


if __name__ == "__main__":
    import asyncio
    from tqdm.auto import tqdm
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")

    # add subcommand 'help'
    help_parser = subparsers.add_parser("help")
    help_parser.add_argument(
        "command", type=str, help="Subcommand to show help for")

    # add subcommand 'download'
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

    # add subcommand 'pdb2fasta'
    pdb2fasta_parser = subparsers.add_parser("pdb2fasta")
    pdb2fasta_parser.add_argument(
        "--fasta_file",
        "-o",
        type=str,
        required=True,
        help="Path to the output fasta file")
    pdb2fasta_parser.add_argument(
        "--pdb_files",
        "-i",
        nargs="+",
        required=True,
        help="Paths to the PDB files to be converted")
    pdb2fasta_parser.add_argument(
        "--multimer_mode",
        "-m",
        choices=["seperate", "joint"],
        default="joint",
        help="Mode of the conversion.")
    pdb2fasta_parser.add_argument(
        "--joint_sep",
        "-s",
        type=str,
        default=":",
        help="Separator of the joint sequence.")

    args = parser.parse_args()

    if args.subcommand == "help":
        subparsers.choices[args.command].print_help()

    elif args.subcommand == "download":
        if args.pdb_ids is not None:
            pdb_ids = args.pdb_ids
        else:
            with open(args.pdb_id_file, "r") as fp:
                pdb_ids = [line.strip() for line in fp]

        for pid in pdb_ids:
            if len(pid) != 4:
                raise ValueError(f"Invalid pdbid {pid}")

        output_path = Path(args.target_path).resolve()
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        elif not output_path.is_dir():
            raise FileNotFoundError(
                f"Target path {args.target_path} is not a directory")

        process_bar = tqdm(total=len(pdb_ids), desc="Downloading PDB files")
        download_failed = []

        def callback(*args):
            if len(args) > 0:
                download_failed.append(args[0])
            process_bar.update(1)

        tasks = [async_download_PDB(pdb_id, str(
            output_path), callback) for pdb_id in pdb_ids]
        loop = asyncio.get_event_loop()
        loop.run_until_complete(asyncio.gather(*tasks))
        process_bar.close()

    elif args.subcommand == "pdb2fasta":
        pdb2fasta(
            args.fasta_file,
            *args.pdb_files,
            multimer_mode=args.multimer_mode,
            joint_sep=args.joint_sep)

    else:
        parser.print_help()
