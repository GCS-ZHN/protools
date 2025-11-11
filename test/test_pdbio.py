import pytest
import requests

from pathlib import Path
from protools import pdbio, utils
from . import tools

def test_get_structure():
    c1 = pdbio.get_structure('data/gen.pdb')[0]['A']
    c2 = pdbio.get_structure('data/4zqk.cif')[0]['B']
    c3 = pdbio.get_structure('data/4ozi.pdb.gz')[0]['H']
    c4 = pdbio.get_structure('data/4zqk.cif.gz')[0]['B']
    assert len(c4) == len(c2)


def test_save_pdb(tmp_path):
    structure = pdbio.get_structure('data/4ozi.pdb.gz')
    out_file = tmp_path / 'out.pdb'
    pdbio.save_pdb(out_file, structure, remarks={
        220: ['Test save_pdb function']
    },
    seqres={
        'A': 'ACDEFGHIKLMNPQRSTVWY',
        'B': 'WYACDEFGHIKLMNPQRSTV'
    })
    structure2 = pdbio.get_structure(out_file)
    assert len(structure2[0]['H']) == len(structure[0]['H'])


@pytest.mark.parametrize(
        'pdb_file, exp_md5sum',
        [
            ('data/4I77.pdb', 'd9b2360744db93ef646c7d4f1b6acf88')
        ]
)
def test_pdb2df(tmp_path: Path, pdb_file: str, exp_md5sum: str):
    m = pdbio.get_structure(pdb_file)
    df = pdbio.pdb2df(m)
    df.to_csv(tmp_path / 'test.csv')
    assert tools.md5_equal(tmp_path / 'test.csv', exp_md5sum)


@pytest.mark.parametrize(
        'pdb_file, exp_md5sum',
        [
            ('data/4I77.pdb', '4304ac6ad0e0d25b019c43304ba30bcd')
        ]
)
def test_read_residue(tmp_path: Path, pdb_file: str, exp_md5sum: str):
    df = pdbio.read_residue(pdb_file)
    df.to_csv(tmp_path / 'test.csv')
    assert tools.md5_equal(tmp_path / 'test.csv', exp_md5sum)


@pytest.mark.parametrize(
        'pdb_id, exp_md5sum',
        [
            ('4i77', '7a9b159c58572f974ed74616b046e5fc'),
            ('4ozi', '3cfca0c847e7d2b45edfc194b18f559a'),
            ('3inj', 'a5c3950eb76e83dd7050f7b8a0ff2762')
        ]
)
@utils.max_retry(err_types=(requests.exceptions.SSLError))
def test_fetch(tmp_path: Path, pdb_id: str, exp_md5sum: str):
    assert tools.md5_equal(pdbio.fetch(pdb_id, tmp_path), exp_md5sum)


@pytest.mark.xfail(reason="Not implemented")
def test_is_aa():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_get_aa_sequence():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_get_aa_residues():
    # TODO
    assert False



@pytest.mark.xfail(reason="Not implemented")
def test_pdb2fasta():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_read_pdb_seq():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_pdb2seq():
    # TODO
    assert False

@pytest.mark.xfail(reason="Not implemented")
def test_read_remarks():
    # TODO
    assert False

@pytest.mark.xfail(reason="Not implemented")
def test_read_missing_residues():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_generate_missing_residues_remarks():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_read_missing_atoms():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_generate_missing_atoms_remarks():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_coord2chain():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_write_modified_residues():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_read_modified_residues():
    # TODO
    assert False
