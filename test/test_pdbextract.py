from pathlib import Path
import pytest

from protools import pdbextract, pdbio
from . import tools


@pytest.mark.parametrize(
        'resi, exp_result',
        [
            ('H1A-4,F6,M8-10', [('H', ('1A', '4')), ('F', ('6', '6')), ('M', ('8', '10'))]),
            ('A1', [('A', ('1', '1'))]),
            ('C2-8B', [('C', ('2', '8B'))])
        ]
)
def test_parse_resi(resi: str, exp_result: list):
    assert list(pdbextract.parse_resi(resi)) == exp_result


@pytest.mark.parametrize(
        'seq, pdb_file, findall, num_find',
        [
            ('QVT', 'data/4I77.pdb', False, 1),
            ('A', 'data/4I77.pdb', True, 33),
            ('A', 'data/4I77.pdb', False, 1),
            ('LL', 'data/4I77.pdb', True, 3),
            ('DIV', 'data/4I77.pdb', True, 1),
        ]
)
def test_find_seq(seq: str, pdb_file: str, findall: bool, num_find: int):
    from Bio.Data import IUPACData
    model = pdbio.get_structure(pdb_file)[0]
    resn_df = pdbextract.find_seq(seq, model, findall=findall)
    resn_df = resn_df.drop_duplicates(
        pdbio.RESIDUE_LEVEL_COLS)
    finded_seq = ''.join([IUPACData.protein_letters_3to1[n.capitalize()] for n in resn_df['resn']])
    mod, div = len(finded_seq) % len(seq), len(finded_seq) // len(seq)
    assert mod == 0
    assert div == num_find
    assert seq * div == finded_seq


@pytest.mark.parametrize(
        'pdb_file',
        ['data/3inj.pdb', 'data/4I77.pdb']
)
def test_get_interface(pdb_file: str):
    pdbextract.get_interface(pdbio.get_structure(pdb_file))


@pytest.mark.parametrize(
        'pdb_file, resi, remain, exp_md5sum',
        [
            ('data/4I77.pdb', 'H1-113,L1-107', True, 'aaf8d8a7b18dbd6d3dcedf256659f4c5')
        ]
)
def test_extract(tmp_path: Path, pdb_file: str, resi: str, remain: True, exp_md5sum: str):
    pdbextract.extract(
        pdb_file, tmp_path / 'extracted.pdb', resi_selection=resi, remain=remain)
    
    assert tools.md5_equal(tmp_path / 'extracted.pdb', exp_md5sum)
