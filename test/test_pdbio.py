from protools import pdbio


def test_get_structure():
    s = pdbio.get_structure('data/gen.pdb')
    s = pdbio.get_structure('data/4zqk.cif')