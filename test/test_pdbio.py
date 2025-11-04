from protools import pdbio


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