#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module for reading and writing PDB files.

PDB standard: https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
"""

from io import StringIO
import warnings
import numpy as np
import requests
import gzip
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union
from collections import OrderedDict
from itertools import product

import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import PDBData, IUPACData

from .seqio import save_fasta, read_seqres, Fasta
from .typedef import FilePathType, FilePathOrIOType, StructureFragmentAAType, StructureFragmentType
from .utils import ensure_path, ensure_fileio
from tqdm.auto import tqdm


__all__ = [
    'is_aa',
    'get_aa_residues',
    'get_aa_sequence',
    'get_structure',
    'save_pdb',
    'download_PDB',
    'async_download_PDB',
    'read_pdb_seq',
    'pdb2fasta',
    'pdb2df']


BACKBONE_ATOMS = ('N', 'CA', 'C', 'O')


def _basic_pdb_column_format(line: str) -> str:
    line = line.strip()
    if len(line) > 80:
        raise ValueError("Line length should not exceed 80 according to the PDB format standard")
    return "{:<80s}\n".format(line)


def is_aa(residue: Residue) -> bool:
    """
    Judge whether a residue is an amino acid.

    Parameters
    ----------
    residue : Residue
        Residue to be judged.

    Returns
    ----------
    is_aa : bool
        Whether the residue is an amino acid.
    """
    if residue.get_resname() in PDBData.protein_letters_3to1_extended:
        if any(map(lambda x: x not in residue, BACKBONE_ATOMS)):
            warnings.warn(
                f"Residue {residue.get_resname()} at {residue.id} not found some backbone atoms.")
        return True
    if all(map(lambda x: x in residue, BACKBONE_ATOMS)):
        warnings.warn(
f"Residue {residue.get_resname()} has all amino acid backbone atom name, \
but is not in the amino acid list. We assume it is an amino acid, \
but you should check the PDB file.")
        return True
    return False


def get_aa_residues(chain: Chain) -> OrderedDict:
    """
    Get the amino acid residues of a chain.
    If multiple residues have the same residue ID,
    the residue with higher occupancy will be used.

    Parameters
    ----------
    chain : Chain
        Chain to be processed.

    Returns
    ----------
    seq : OrderedDict
        Amino acid residues of the chain.
    """
    seq = OrderedDict()
    for residue in chain:
        if is_aa(residue):
            resi = residue.get_id()[1:]
            if resi in seq:
                warnings.warn(
f"Chain {chain.get_id()} at residue {resi} has multiple amino acid residues: \
{residue.get_resname()} and {seq[resi].get_resname()}. The residue with higher \
occupancy will be used.")
                prev_occ = [atom.get_occupancy() for atom in seq[resi]]
                curr_occ = [atom.get_occupancy() for atom in residue]
                prev_occ = sum(prev_occ) / len(prev_occ)
                curr_occ = sum(curr_occ) / len(curr_occ)
                if curr_occ > prev_occ:
                    seq[resi] = residue
            else:
                seq[resi] = residue
        else:
            continue
    return seq


def get_aa_sequence(chain: Chain, standard: bool = True, unknown_aa: str = 'X') -> Seq:
    """
    Get the amino acid sequence of a chain.

    Parameters
    ----------
    chain : Chain
        Chain to be processed.
    standard : bool, optional
        Whether only standard amino acids are included, by default True.
    unknown_aa : str, optional
        Symbol to represent unknown amino acids, by default 'X'.

    Returns
    ----------
    seq : Seq
        Amino acid sequence of the chain.
    """
    aa_3to1 = PDBData.protein_letters_3to1_extended
    if standard:
        aa_3to1 = PDBData.protein_letters_3to1
    seq = get_aa_residues(chain)
    return Seq(''.join(aa_3to1.get(residue.get_resname(), unknown_aa) for residue in seq.values()))


def get_structure(pdb_file: FilePathType, structure_id: str = 'pdb') -> Structure:
    """
    Get the structure of a PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.
    structure_id : str, optional
        ID of the structure, by default 'pdb'.

    Returns
    ----------
    structure : Structure
        Structure of the PDB file.
    """
    pdb_file = ensure_path(pdb_file)
    if not pdb_file.exists():
        raise FileNotFoundError(f"Could not find PDB file {pdb_file}")

    if pdb_file.suffix.lower() == '.cif':
        parser = MMCIFParser(QUIET=True)
    elif pdb_file.suffix.lower() == '.pdb':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file type {pdb_file.suffix}")
    return parser.get_structure(structure_id, pdb_file)


def _write_seqres(target_path: FilePathOrIOType, seqres: dict):
    """
    Write SEQRES records to a PDB file.

    Parameters
    ----------
    target_path : str
        Path to the output PDB file.
    seqres : dict
        SEQRES records to be written to the PDB file.
        The key is the chain ID and the value is the
        one-letter amino acid sequence or a list of
        three-letter amino acid sequence or nucleotide

    Notes
    ----------
    SEQRES records maybe amino acid or nucleotide sequence.
    But this function only supports amino acid sequence if
    the sequence is a one-letter amino acid string. If you
    want to write nucleotide sequence, you should input a list
    of three-letter nucleotide.
    """
    format_line = lambda k: 'SEQRES {:>3d} {:1s} {:>4d} ' + ''.join([' {:>3s}'] * k)
    f, need_close = ensure_fileio(target_path, 'w')
    for chain, seq in seqres.items():
        if isinstance(seq, (str, Seq)):
            seq = [PDBData.protein_letters_1to3.get(aa, 'UNK') for aa in seq]
        elif isinstance(seq, (list, tuple)):
            pass
        else:
            raise ValueError("Invalid SEQRES format")
        for i in range(0, len(seq), 13):
            seq_slice = seq[i:i+13]
            line = format_line(len(seq_slice)).format(i//13+1, chain, len(seq), *seq_slice)
            f.write(_basic_pdb_column_format(line))
    if need_close:
        f.close()
    else:
        f.flush()


def _write_remark(target_path: FilePathOrIOType, remark_id: int, *remarks: str):
    format_line = lambda k: 'REMARK {:>3d} {:s}'
    f, need_close = ensure_fileio(target_path, 'w')
    if remark_id < 0 or remark_id > 999:
        raise ValueError("Remark ID should be between 0 and 999")
    for line in remarks:
        f.write(_basic_pdb_column_format(format_line(remark_id).format(remark_id, line)))
    if need_close:
        f.close()
    else:
        f.flush()


def read_modified_residues(pdbfile: Path) -> pd.DataFrame:
    """
    Read modified residues from MODRES Record
    """
    result = {
        'id_code': [],            # column 8-11
        'res_name': [],           # column 13-15
        'chain_id': [],           # column 17
        'sequence_number' : [],   # column 19-22
        'insertion_code': [],     # column 23
        'standard_res_name': [],  # column 25-27
        'comment': []             # column 30-
    }
    with open(pdbfile) as f:
        for line in f:
            if line.startswith('MODRES'):
                id_code = line[7:11].strip()
                res_name = line[12:15].strip()
                chain_id = line[16]
                ssseq = int(line[18:22].strip())
                inscode = line[22]
                standard_res_name = line[24:27].strip()
                comment = line[29:].strip()
                result['id_code'].append(id_code)
                result['res_name'].append(res_name)
                result['chain_id'].append(chain_id)
                result['sequence_number'].append(ssseq)
                result['insertion_code'].append(inscode)
                result['standard_res_name'].append(standard_res_name)
                result['comment'].append(comment)
    return pd.DataFrame(result)


def write_modified_residues(data: pd.DataFrame, target_path: FilePathOrIOType):
    """
    Generate MODRES Record
    """
    if len(data) == 0:
        return
    format_line = 'MODRES {:<4s} {:<3s} {:<1s} {:>4d}{:<1s} {:<3s}  {:<s}'
    f, need_close = ensure_fileio(target_path, 'w')
    for _, row in data.iterrows():
        f.write(_basic_pdb_column_format(format_line.format(*row.values)))
    if need_close:
        f.close()
    else:
        f.flush()


def save_pdb(
        output_path: FilePathType, 
        *entities: StructureFragmentAAType, 
        remarks: Optional[Dict[int, Iterable[str]]] = None,
        seqres: dict = None,
        modres: pd.DataFrame = None) -> None:
    """
    Save entities to a PDB file.

    Parameters
    ----------
    output_path : str
        Path to the output PDB file.
    entities : Structure, Model, Chain, or Residue
        Entities to save.
    remarks : Dict[int, Iterable[str]], optional
        Remarks to be written to the PDB file.
    seqres : dict, optional
        SEQRES records to be written to the PDB file.
        The key is the chain ID and the value is the
        one-letter amino acid sequence. Only works
        when the entities are not Residue objects.
    modres : pd.DataFrame, optional
        MODRES records to be written to the PDB file.
        The DataFrame should have the following columns:
        id_code, res_name, chain_id, sequence_number,
        insertion_code, standard_res_name, comment.

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
    
    output_path = ensure_path(output_path)

    def _loop_type_check(entities, *allowed_types, allow_none=False):
        for entity in entities:
            if not isinstance(entity, allowed_types):
                if not allow_none or entity is not None:
                    raise ValueError(f"Unsupported type {type(entity)}")
            yield entity

    if isinstance(entities[0], Structure):
        if len(entities) > 1:
            raise RuntimeError("Only accept one Structure object")
        pdb_io = PDBIO()
        pdb_io.set_structure(entities[0])
        with open(output_path, "w") as fp:
            _write_remark(fp, 220, "Created by protools")
            if remarks is not None:
                for k, v in remarks.items():
                    _write_remark(fp, k, *v)
            if seqres is not None:
                _write_seqres(fp, seqres)
            if modres is not None:
                write_modified_residues(modres, fp)
            pdb_io.save(fp)

    elif isinstance(entities[0], Model):
        structure = Structure("pdb")
        if len(entities) > 1:
            # TODO: support multiple models with differnt SEQRES
            raise RuntimeError("Multiple Model objects are not implemented yet")
        for model in _loop_type_check(entities, Model):
            structure.add(model)
        save_pdb(
            output_path,
            structure,
            remarks=remarks,
            seqres=seqres,
            modres=modres)

    elif isinstance(entities[0], Chain):
        models = [Model("model_0")]
        for chain in _loop_type_check(entities, Chain, allow_none=True):
            if chain is None:
                models.append(Model(f"model_{len(models)}"))
            else:
                models[-1].add(chain)
        save_pdb(
            output_path,
            *models,
            remarks=remarks,
            seqres=seqres,
            modres=modres)

    elif isinstance(entities[0], Residue):
        if seqres is not None:
            raise ValueError("SEQRES is not supported for Residue object")
        chains = [Chain("A")]
        for residue in _loop_type_check(entities, Residue, allow_none=True):
            if residue is None:
                chain_id = chr(ord(chains[-1].get_id()) + 1)
                if chain_id > "Z":
                    raise ValueError("Too many chains")
                chains.append(Chain(chain_id))
            else:
                chains[-1].add(residue)
        save_pdb(
            output_path,
            *chains,
            remarks=remarks,
            modres=modres)

    else:
        raise TypeError(f"Unsupported type {type(entities[0])}")


def fetch(pdb_id: str, target_dir: FilePathType, server: str = 'https://files.rcsb.org') -> Path:
    """
    Download a PDB file from the PDB database.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the PDB file to download.

    target_dir : str
        Path to save the downloaded PDB file.

    server : str, optional
        URL of the PDB server to download from.

    Returns
    ----------
    file : Path
        Path to the downloaded PDB file.

    Raises
    ----------
    FileNotFoundError
        If the PDB file could not be downloaded.
    """
    target_dir = ensure_path(target_dir)
    target_dir.mkdir(exist_ok=True, parents=True)
    file_name_format = {
        'pdb': f'{pdb_id[1:3]}/pdb{pdb_id}.ent.gz',
        'mmCIF': f'{pdb_id[1:3]}/{pdb_id}.cif.gz'
    }
    type2suffix = {
        'pdb': '.pdb',
        'mmCIF': '.cif'
    }
    tmp_file = target_dir / f'{pdb_id}.tmp'

    for pdb_dir, pdb_type in product(['divided', 'obsolete'], type2suffix.keys()):
        try:
            url = f'{server}/pub/pdb/data/structures/{pdb_dir}/{pdb_type}/{file_name_format[pdb_type]}'
            response = requests.get(url, stream=True)
            response.raise_for_status()
            total = int(response.headers.get('content-length', 0))
            with tqdm.wrapattr(response.raw, 'read', total=total,
                               desc=f'Downloading {pdb_id}') as stream:
                with gzip.open(stream, mode='r') as f:
                    with tmp_file.open('wb') as out:
                        for chunk in iter(lambda: f.read(1024), b''):
                            out.write(chunk)
            path = target_dir / (pdb_id + type2suffix[pdb_type])
            tmp_file.rename(path)
            return path
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                print(f'Failed to fetch {pdb_id} [{pdb_dir}/{pdb_type}]: {e}, try other')
                continue
            raise e

    raise FileNotFoundError(f'Failed to fetch {pdb_id}')


def read_pdb_seq(entity: StructureFragmentType, standard: bool = True, unknown_aa:str = 'X') -> Iterable[Tuple[str, str, Seq]]:
    """
    Extract the sequence of a S.

    Parameters
    ----------
    entity : StructureFragmentType
        Structure, Model or Chain object.
    standard : bool, optional
        Whether only standard amino acids are included, by default True.
    unknown_aa : str, optional
        Symbol to represent unknown amino acids, by default 'X'.

    Returns
    ----------
    seq_iter : Iterable[Tuple[str, str, Seq]]
        An iterable of tuples (model_id, chain_id, seq).
    """
    if entity.level not in ['S', 'M', 'C']:
        raise ValueError("Require Structure, Model or Chain object")
    
    if entity.level == 'C':
        chains = [entity]
    else:
        chains = entity.get_chains()

    for chain in chains:
        model = chain.get_parent()
        model_id = model.get_id() if model else 0
        seq = get_aa_sequence(chain, standard=standard, unknown_aa=unknown_aa)
        yield model_id, chain.get_id(), seq


def pdb2fasta(
        fasta_path: str,
        *pdb_files: str,
        multimer_mode: str = 'joint',
        selected_chains: Optional[str] = None,
        standard: bool = True,
        unknown_aa: str = 'X',
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

    standard : bool, optional
        Whether only standard amino acids are included, by default True.

    unknown_aa : str, optional
        Symbol to represent unknown amino acids, by default 'X'.

    selected_chains: str, optional
        If provided, only selected chain will be
        included.

    joint_sep : str, optional
        Separator of the joint sequence. The default is ':'.
        Only works when multimer_mode is 'joint'.
        Different tools may have different requirements
        to process the joint sequence.

    """
    if len(pdb_files) == 0:
        raise ValueError("No PDB files to convert")

    fasta_path = ensure_path(fasta_path)

    def _iter():
        for pdb_file in pdb_files:
            pdb_id = Path(pdb_file).stem
            seq_iter = read_pdb_seq(get_structure(pdb_file), standard=standard, unknown_aa=unknown_aa)
            seq_iter = filter(lambda x: selected_chains is None or x[1] in selected_chains, seq_iter)
            if multimer_mode == 'seperate':
                for model_id, chain_id, seq in seq_iter:
                    seq_id = f"{pdb_id}_{model_id}_{chain_id}"
                    yield SeqRecord(seq, id=seq_id, description='')

            elif multimer_mode == 'joint':
                seq = joint_sep.join(str(seq) for _, _, seq in seq_iter)
                seq_id = pdb_id
                yield SeqRecord(Seq(seq), id=seq_id, description='')

            else:
                raise ValueError(
                    f"Unsupported multimer mode {multimer_mode}")
    save_fasta(_iter(), fasta_path, mkdir=True)


def pdb2seq(path: FilePathType) -> Fasta:
    """
    Convert a PDB file to a Fasta object.

    Parameters
    ----------
    path : str
        Path to the PDB file.

    Returns
    ----------
    fasta : Fasta
        Fasta object of the PDB file.
    """
    seqs = read_seqres(path)
    if len(seqs) == 0:
        seqs = read_pdb_seq(get_structure(path))
        seqs = Fasta((x[1], x[2]) for x in seqs)
    return seqs


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
            'seqid': resids[1],
            'inscode': resids[2],
            'chain': chain.get_id(), 
            'x': coord[0], 
            'y': coord[1], 
            'z': coord[2],
            'element': atom.element,
            'model': model.get_id()
        }
        for attr in extra_attrs:
            try:
                res[attr] = getattr(atom, attr)
            except AttributeError as e:
                try:
                    res[attr] = getattr(atom, f'get_{attr}')()
                except:
                    raise e
        return res
    
    return pd.DataFrame(_atom_to_dict(atom) for atom in entity.get_atoms())


def read_residue(pdb: Union[FilePathType, str, Entity], mode='centroid') -> pd.DataFrame:
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
    model chain seqid resn                                
    0     A     100  GLY   7.494460 -33.223431  16.242475
                101  ARG   3.844297 -34.536686  13.638067
                102  ASP   2.507354 -39.138135  16.654611
                103  THR   3.519204 -40.032841  12.495594
                104  SER   1.635237 -43.545750  14.313151
    """
    if isinstance(pdb, Entity):
        structure = pdb
    elif isinstance(pdb, (str, Path)):
        structure = get_structure(pdb)
    else:
        raise TypeError(f"Unsupported type {type(pdb)}")

    if mode == 'centroid':
        df = pdb2df(structure, 'mass')
        df['x'] *= df['mass']
        df['y'] *= df['mass']
        df['z'] *= df['mass']
        df = df.drop(['id', 'name', 'element'], axis=1)
        df = df.groupby(['model', 'chain', 'seqid', 'resn']).sum()
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
        df = df.groupby(['model', 'chain', 'seqid', 'resn']).mean()
        return df
    
    raise ValueError('mode must be one of "centroid", "fuc", "CA"')


def read_remarks(pdbfile: Path, ignore_non_standard: bool = False) -> Dict[int, str]:
    """
    Read all remarks defined in PDB.

    Parameters
    ----------
    pdbfile : Path
        Path to the PDB file.
    ignore_non_standard : bool, optional
        If True, ignore remarks that are not standard PDB remarks.

    Returns
    ----------
    remarks : Dict[int, str]
        A dictionary where the keys are the remark type codes
    """
    remarks = {}
    with open(pdbfile) as f:
        for line in f:
            if line.startswith('REMARK'):
                if ignore_non_standard and not line[7:10].isdigit():
                    continue
                remark_type_code = int(line[7:10])
                remark_content = line[11:]
                if remark_type_code not in remarks:
                    remarks[remark_type_code] = remark_content
                else:
                    remarks[remark_type_code] += remark_content
    return remarks


def read_missing_residues(pdbfile: Path) -> pd.DataFrame:
    """
    Read missing residues from REMARK 465
    """
    remark = read_remarks(pdbfile).get(465, "")
    remark = StringIO(remark)
    reach_title = False
    result = {
        'model': [],
        'res_name': [],
        'chain_id': [],
        'sequence_number' : [],
        'insertion_code': [],
    }
    for line in remark:
        if line.strip() == 'M RES C SSSEQI':
            reach_title = True
            continue
        if not reach_title:
            continue
        model = line[2]
        model = int(model) if model.isdigit() else 0
        resn = line[4:7].strip()
        chain_id = line[8]
        ssseq = int(line[10:15].strip())
        inscode = line[15]
        result['model'].append(model)
        result['res_name'].append(resn)
        result['chain_id'].append(chain_id)
        result['sequence_number'].append(ssseq)
        result['insertion_code'].append(inscode)
    return pd.DataFrame(result)


def generate_missing_residues_remarks(data: pd.DataFrame) -> Dict[int, list]:
    """
    Generate REMARK 465
    """
    remarks = [
        '',
        'MISSING RESIDUES',
        'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE',
        'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN',
        'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)',
        '',
        '  M RES C SSSEQI'
    ]
    format_line = '  {:>1} {:>3s} {:>1s} {:>5d}{:>1s}'
    for _, row in data.iterrows():
        remarks.append(format_line.format(
            row['model'],
            row['res_name'],
            row['chain_id'],
            row['sequence_number'],
            row['insertion_code'],
        ))
    return {465: remarks}


def read_missing_atoms(pdbfile: Path) -> pd.DataFrame:
    """
    Read missing atoms from REMARK 470
    """
    remark = read_remarks(pdbfile).get(470, "")
    remark = StringIO(remark)
    reach_title = False
    result = {
        'model': [],
        'res_name': [],
        'chain_id': [],
        'sequence_number' : [],
        'insertion_code': [],
        'atom_name': [],
    }
    for line in remark:
        if line.strip() == 'M RES CSSEQI  ATOMS':
            reach_title = True
            continue
        if not reach_title:
            continue
        model = line[2]
        model = int(model) if model.isdigit() else 0
        resn = line[4:7].strip()
        chain_id = line[8]
        ssseq = int(line[9:13].strip())
        inscode = line[13]
        atom_names = line[16:].strip().split()
        result['model'].append(model)
        result['res_name'].append(resn)
        result['chain_id'].append(chain_id)
        result['sequence_number'].append(ssseq)
        result['insertion_code'].append(inscode)
        result['atom_name'].append(atom_names)
    return pd.DataFrame(result)


def generate_missing_atoms_remarks(data: pd.DataFrame) -> Dict[int, list]:
    """
    Generate remark 470
    """
    remarks = [
        '',
        'MISSING ATOMS',
        'THE FOLLOWING RESIDUES HAVE MISSING ATOMS(M=MODEL NUMBER;',
        'RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;',
        'I=INSERTION CODE):',
        '  M RES CSSEQI  ATOMS',
    ]
    format_func = lambda k: '  {:>1} {:>3s} {:>1s}{:>4d}{:>1s} ' + ''.join(['  {:<3s}'] * k)
    for _, row in data.iterrows():
        format_line = format_func(len(row['atom_name']))
        remarks.append(format_line.format(
            row['model'],
            row['res_name'],
            row['chain_id'],
            row['sequence_number'],
            row['insertion_code'],
            *row['atom_name'],
        ))
    return {470: remarks}


def coord2chain(
        coord: np.ndarray, seq: str, chain_id: str = 'A', atoms: List[str] = None,
        init_residue_number: int = 1, init_serial_number: int = 1) -> Chain:
    """
    Convert coordinates to a Chain object. It's useful to build a new structure
    from scratch.

    Parameters
    ----------
    coord : np.ndarray
        Coordinates of the atoms.
    seq : str
        Sequence of the chain.
    chain_id : str, optional
        ID of the chain, by default 'A'.
    atoms : List[str], optional
        List of atom names, by default None.
    init_residue_number : int, optional
       Initial residue number, by default 1.
    init_serial_number : int, optional
        Initial serial number, by default 1.


    Returns
    ----------
    chain : Chain
        Chain object with the coordinates and sequence.
    """
    if atoms is None:
        atoms = BACKBONE_ATOMS
    chain = Chain(chain_id)
    coord = coord.reshape(len(seq), len(atoms), 3).astype(np.float32)
    serial_number = init_serial_number
    for res_idx, res_letter in enumerate(seq):
        res_name = IUPACData.protein_letters_1to3[res_letter.upper()].upper()
        residue = Residue((' ', res_idx + init_residue_number, ' '), res_name, '    ')
        chain.add(residue)
        residue_coord = coord[res_idx]
        for atom_idx, atom_name in enumerate(atoms):
            atom = Atom(
                name=atom_name,
                coord=residue_coord[atom_idx],
                bfactor=0.0,
                occupancy=1.0,
                altloc=' ',
                fullname=f' {atom_name:<3s}',
                serial_number=serial_number,
                element=atom_name[0],
            )
            residue.add(atom)
            serial_number += 1
    return chain


if __name__ == "__main__":
    import asyncio
    from argparse import ArgumentParser

    from tqdm.auto import tqdm
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    # add subcommand 'help'
    help_parser = subparsers.add_parser("help")
    help_parser.add_argument(
        "command", type=str, help="Subcommand to show help for")

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
        "--chain-ids",
        "-c",
        help="Selected chains will be extracted."
    )
    pdb2fasta_parser.add_argument(
        "--joint_sep",
        "-s",
        type=str,
        default=":",
        help="Separator of the joint sequence.")

    args = parser.parse_args()

    if args.cmd == "help":
        subparsers.choices[args.command].print_help()

    elif args.cmd == "pdb2fasta":
        pdb2fasta(
            args.fasta_file,
            *args.pdb_files,
            multimer_mode=args.multimer_mode,
            selected_chains=args.chain_ids,
            joint_sep=args.joint_sep)

    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")
