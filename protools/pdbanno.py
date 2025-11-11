import logging
import re
import asyncio
import subprocess
from itertools import product, starmap
from multiprocessing import Pool
from pathlib import Path
from typing import Union, Tuple, Iterable

import pandas as pd
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.Data import IUPACData

from protools import pdbio
from protools.aa import RESIDUE2SIDECHAIN_3LETTER
from protools.typedef import FilePathType, StructureFragmentAAType, StructureFragmentType
from protools.utils import ensure_path, CmdWrapperBase

from scipy.spatial.distance import cdist

_LOGGER = logging.getLogger(__name__)


def neighbor_water_count(
        entity: StructureFragmentType, 
        level: str = 'R', 
        threshold: float = 2.8) -> Union[int, pd.Series]:
    """Count the number of water molecules in a PDB file.

    Parameters
    ----------
    entity : Union[Chain, Model, Structure]
        PDB Entity.
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
    >>> from protools import pdbanno, pdbio
    >>> entity = pdbio.get_structure('3inj.pdb')
    >>> pdbanno.neighbor_water_count(entity)
    model  chain  seqid  resn  inscode
    0      A      7      ALA              0
                8      VAL              0
                9      PRO              0
                10     ALA              0
                11     PRO              1
                                        ..
        H      496    PRO              0
                497    GLN              0
                498    LYS              2
                499    ASN              0
                500    SER              3
    Length: 3953, dtype: int64
    """
    assert isinstance(entity, StructureFragmentType), \
        f"pdb must be a StructureFragmentType"
    df = pdbio.pdb2df(entity)

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
        level_columns = pdbio.RESIDUE_LEVEL_COLS
    else:
        raise ValueError(f"Unknown level: {level}")


    dist = cdist(res_hatom_coord, water_oxygen_coord, 'euclidean') <= threshold
    dist = pd.DataFrame(index=res_hatom_coord.index, data=dist)
    if level == 'S':
        return dist.any().sum()

    dist = dist.join(res_hatom[level_columns])
    return dist.groupby(level_columns).any().sum(axis=1)


def calc_sasa(
        entity: Union[StructureFragmentType, Residue], 
        radius: float = 1.4, 
        standard: bool = True) -> pd.Series:
    """
    calculate the solvent accessible surface area (SASA) of each residue
    in a structure fragment.

    Parameters
    ----------
    entity : StructureFragmentType or Residue
        A Bio.PDB.Structure.Structure, Bio.PDB.Model.Model,
        Bio.PDB.Chain.Chain or Bio.PDB.Residue.Residue object.
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

    Examples
    ----------
    >>> from protools import pdbanno
    >>> pdb = pdbio.get_structure('3inj.pdb')
    >>> pdbanno.calc_sasa(pdb[0])
    model  chain  seqid  resn  inscode
    0      A      7     ALA             139.496754
                  8     VAL              24.683581
                  9     PRO              49.398717
                  10    ALA              63.504719
                  11    PRO               4.285836
                                        ...    
           H      496   PRO              38.512169
                  497   GLN               0.000000
                  498   LYS               0.000000
                  499   ASN               0.000000
                  500   SER               0.000000
    Length: 3953, dtype: float64
    >>> pdbanno.calc_sasa(pdb[0]['H'])
    model  chain  seqid  resn  inscode
    0      H      6     GLN             142.226811
                  7     ALA              64.825545
                  8     VAL              35.011793
                  9     PRO              30.076666
                  10    ALA              72.322366
                                         ...    
                  496   PRO             101.058627
                  497   GLN              52.869222
                  498   LYS              40.743604
                  499   ASN              93.610187
                  500   SER             129.328349
    Length: 495, dtype: float64
    >>> pdbanno.calc_sasa(pdb[0]['H'][6])
    model  chain  seqid  resn  inscode
    0      H      6      GLN              285.154046
    dtype: float64

    Notes
    ----------
    the sasa of a residue is dependent on the full input
    fragment, not just the residue itself. So, if you 
    provide a residue, a chain or a model, the sasa of
    the same residue may be different.
    """
    if entity.level == 'R':
        residues = [entity]
    elif entity.level in ('C', 'M', 'S'):
        residues = list(entity.get_residues())
    else:
        raise ValueError(f"Unknown entity level: {entity.level}")

    sasa_calc = ShrakeRupley(probe_radius=radius)
    sasa_calc.compute(entity, level='R')
    result = dict()

    for res in residues:
        if standard and res.resname.capitalize() not in IUPACData.protein_letters_3to1:
            continue

        if res.resname == 'HOH':
            continue
        
        chain = res.get_parent()
        chain_id = chain.get_id() if chain else 'A'
        model = chain.get_parent() if chain else None
        model_id = model.get_id() if model else 0
        resid = res.get_id()[1]
        inscode = res.get_id()[2]
        index = (model_id, chain_id, resid, res.resname, inscode)
        result[index] = res.sasa
    res = pd.Series(result)
    res.index.set_names(pdbio.RESIDUE_LEVEL_COLS, inplace=True)
    return res


def _calc_sasa_from_pdb(
        pdb_file: FilePathType,
        output_dir: FilePathType,
        model_idx: int = 0):
    """
    A wrapper for calc_sasa. just
    used for `calc_sasa_from_pdbs`.

    Parameters
    ----------
    pdb_file : FilePathType
        Path to the PDB file.
    output_dir : FilePathType
        Path to the output directory.
    model_idx : int, optional
        The index of the model to be used, by default 0.
    """
    _LOGGER.info(f"calculating sasa for {pdb_file}")
    pdb_file = ensure_path(pdb_file)
    output_dir = ensure_path(output_dir)
    structure = pdbio.get_structure(pdb_file)
    sasa = calc_sasa(structure[model_idx])
    output_file = output_dir / f"{pdb_file.stem}.csv"
    sasa.to_csv(output_file)
    _LOGGER.info(f"done: {output_file}")


def calc_sasa_from_pdbs(
        pdb_dir: FilePathType, 
        output_dir: FilePathType, 
        model_idx: int = 0, 
        num_worker: int = 1):
    """
    SASA batch calculation for PDB files.

    Parameters
    ----------
    pdb_dir : FilePathType
        Path to the directory containing PDB files.
    output_dir : FilePathType
        Path to the output directory.
    model_idx : int, optional
        The index of the model to be used, by default 0.
    num_worker : int, optional
        Number of workers, by default 1.
    """
    pdb_dir = ensure_path(pdb_dir)
    output_dir = ensure_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    schedules = product(pdb_dir.glob('*.pdb'), [output_dir], [model_idx])

    if num_worker == 1:
        list(starmap(_calc_sasa_from_pdb, schedules))

    elif num_worker > 1:
        with Pool(num_worker) as p:
            p.starmap(_calc_sasa_from_pdb, schedules)

    else:
        raise ValueError("num_worker must be greater than 0")
    
    _LOGGER.info("all done")


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
    entity1_df = pdbio.pdb2df(entity1)
    entity2_df = pdbio.pdb2df(entity2)
    entity1_df.set_index(pdbio.RESIDUE_LEVEL_COLS, inplace=True)
    entity2_df.set_index(pdbio.RESIDUE_LEVEL_COLS, inplace=True)

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

class TMalign(CmdWrapperBase):
    """
    TMalign python binding.
    """
    def __init__(self, 
                 num_workers: int = 1,
                 byresi: int = 0,
                 clean: bool = True,
                 tmalign_path: str = 'TMalign'):
        super().__init__(
            tmalign_path,
            short_mode=True,
            num_workers=num_workers,
            install_help="Please install TMalign from https://zhanggroup.org/TM-align/")
        self.byresi = byresi
        self.clean = clean

    def _preprocess_(self,
                    query_pdb: FilePathType,
                    native_pdb: FilePathType,
                    output_dir: FilePathType, prefix: str = None):
        query_pdb = ensure_path(query_pdb)
        native_pdb = ensure_path(native_pdb)
        output_dir = ensure_path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if prefix is None:
            prefix = query_pdb.stem
        else:
            if prefix != Path(prefix).stem:
                raise ValueError(f"Invalid prefix: {prefix}")
        return query_pdb, native_pdb, output_dir, prefix
    
    def _postprocess_(self, 
                     result: subprocess.CompletedProcess, 
                     output_dir: Path, prefix: str = None):
        result = result.stdout.decode('utf-8')
        tm_score = float(re.search(r"TM-score=\s+([\d\.]+)", result).group(1))
        rmsd = float(re.search(r"RMSD=\s+([\d\.]+)", result).group(1))

        if self.clean:
            for suffix in ['', '_all', '_all_atm', '_atm', '_all_atm_lig']:
                file = output_dir / (prefix + suffix)
                file_pml = output_dir / (prefix + suffix + '.pml')
                if file.exists():
                    file.unlink()
                if file_pml.exists():
                    file_pml.unlink()
        return tm_score, rmsd

    def __call__(self, 
                 query_pdb: FilePathType,
                 native_pdb: FilePathType,
                 output_dir: FilePathType,
                 prefix: str = None) -> Tuple[float, float]:
        """
        Align two pdb files using TMalign and return TM-score.

        Parameters
        ----------
        query_pdb : FilePathType
            Path to query pdb file.

        native_pdb : FilePathType
            Path to native pdb file.

        output_dir : FilePathType
            Path to output directory.
        prefix : str
            Prefix of output file.

        Returns
        ----------
        Tuple[float]
            TM-score and RMSD.
        """
        query_pdb, native_pdb, output_dir, prefix = self._preprocess_(
            query_pdb, native_pdb, output_dir, prefix)
        result = super().__call__(
            query_pdb,
            native_pdb,
            o=output_dir / prefix,
            byresi=self.byresi,
            stdout=subprocess.PIPE)
        return self._postprocess_(result, output_dir, prefix)

    async def async_call(self, 
                         query_pdb: FilePathType,
                         native_pdb: FilePathType,
                         output_dir: FilePathType,
                         prefix: str = None) -> Tuple[float, float]:
        """
        Align two pdb files using TMalign and return TM-score.

        Parameters
        ----------
        query_pdb : FilePathType
            Path to query pdb file.

        native_pdb : FilePathType
            Path to native pdb file.

        output_dir : FilePathType
            Path to output directory.
        prefix : str
            Prefix of output file.

        Returns
        ----------
        Tuple[float]
            TM-score and RMSD.
        """
        query_pdb, native_pdb, output_dir, prefix = self._preprocess_(
            query_pdb, native_pdb, output_dir, prefix)
        result = await super().async_call(
            query_pdb,
            native_pdb,
            o=output_dir / prefix,
            byresi=self.byresi,
            stdout=subprocess.PIPE)
        return self._postprocess_(result, output_dir, prefix)

    async def async_batch_align(self,
                    pair_iterable: Iterable[Tuple[str, str]],
                    output_dir: FilePathType) -> pd.DataFrame:
        """
        Align a batch of pdb files using TMalign and return TM-score.

        Parameters
        ----------
        pair_iterable : Iterable[Tuple[str, str]]
            Iterable of pdb file pairs.
        output_dir : FilePathType
            Path to output directory.

        Returns
        ----------
        pd.DataFrame
            TM-score and RMSD for each input pdb file pair.
        """
        output_dir = ensure_path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb_info = list(pair_iterable)
        pdb_info = pd.DataFrame(pdb_info, columns=['gen', 'native'])
        params = map(
            lambda x: (x.gen, x.native, output_dir),
            pdb_info.itertuples(index=False))
        if self.semaphore._value > 1:
            tasks = list(starmap(self.async_call, params))
            results = await asyncio.gather(*tasks)
            results = pd.DataFrame(
                results,
                columns=['tmscore', 'rmsd'])

        else:
            results = pd.DataFrame(
                starmap(self, params),
                columns=['tmscore', 'rmsd'])
        return pd.concat([pdb_info, results], axis=1)
    
    def batch_align(self,
                    pair_iterable: Iterable[Tuple[str, str]],
                    output_dir: FilePathType) -> pd.DataFrame:
        """
        Align a batch of pdb files using TMalign and return TM-score.

        Parameters
        ----------
        pair_iterable : Iterable[Tuple[str, str]]
            Iterable of pdb file pairs.
        output_dir : FilePathType
            Path to output directory.

        Returns
        ----------
        pd.DataFrame
            TM-score and RMSD for each input pdb file pair.

        See Also
        ----------
        async_batch_align

        Notes
        -----------
        This method is not recommended for use in a Jupyter notebook,
        Otherwise, a RuntimeError will be raised. Or you can use
        `nest_asyncio` to avoid this error.
        """
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
        if loop.is_running():
            raise RuntimeError("Event loop already exist and is running, \
please use `async_batch_align` if you want to use asyncio.")
        return loop.run_until_complete(
            self.async_batch_align(pair_iterable, output_dir))


if __name__ == '__main__':
    from argparse import ArgumentParser

    common_parser = ArgumentParser(add_help=False)
    common_parser.add_argument("--pdbfile", "-i", type=Path, required=True)
    common_parser.add_argument("--output", "-o", type=Path, required=True)
    common_parser.add_argument("--model_idx", type=int, default=0)
    common_parser.add_argument("--log_level", "-l", type=str, default="INFO")

    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd')

    sasa_parser = subparsers.add_parser('sasa', parents=[common_parser])
    sasa_parser.add_argument("--num_worker", "-n", type=int, default=1)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.getLevelName(args.log_level),
        format='%(process)d-%(asctime)s-%(levelname)s-%(message)s')
    
    if args.cmd == 'sasa':
        calc_sasa_from_pdbs(
            args.pdbfile, 
            args.output, 
            args.model_idx, 
            args.num_worker)
        
    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")
