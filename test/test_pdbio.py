from protools import pdbio


def test_get_structure():
    c1 = pdbio.get_structure('data/gen.pdb')[0]['A']
    c2 = pdbio.get_structure('data/4zqk.cif')[0]['B']
    c3 = pdbio.get_structure('data/4ozi.pdb.gz')[0]['H']
    c4 = pdbio.get_structure('data/4zqk.cif.gz')[0]['B']
    assert len(c4) == len(c2)
    