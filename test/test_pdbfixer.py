import pytest

from protools import pdbfixer
from .tools import md5_equal


@pytest.mark.parametrize(
        "input_pdb, exp_md5sum, restart_by_chain", 
        [
            ('data/4ozi.pdb', '240b96faea97991fd0bcf747c0e7fe0a', False),
            ('data/4ozi.pdb', '612b2ecd5f1430a60c3987f110ea8378', True)
        ]
)
def test_renumber_residue(tmp_path, input_pdb, exp_md5sum, restart_by_chain):
    """
    Test the renumber_residue function.
    """
    output_pdb = tmp_path / 'renumbered.pdb'
    pdbfixer.renumber_residue(
        pdb_file=input_pdb,
        out_file=output_pdb,
        model_idx=0,
        start=1,
        restart_by_chain=restart_by_chain
    )
    assert md5_equal(output_pdb, exp_md5sum), "Renumbered PDB file does not match expected MD5 checksum"


@pytest.mark.parametrize(
        "input_pdb, exp_md5sum, chain_id_map, chain_order", 
        [
            ('data/4I77.pdb', 'c981b42afad847650a50ebd3ebc06887',{'Z': 'T'}, ['H', 'L', 'Z']),
        ]
)
def test_renumber_residue_chain_id_map(tmp_path, input_pdb, exp_md5sum, chain_id_map, chain_order):
    """
    Test the renumber_residue function with chain_id_map option.
    """
    output_pdb = tmp_path / 'renumbered.pdb'
    pdbfixer.renumber_residue(
        pdb_file=input_pdb,
        out_file=output_pdb,
        model_idx=0,
        start=1,
        chain_id_map=chain_id_map,
        chain_order=chain_order
    )
    assert md5_equal(output_pdb, exp_md5sum), \
        "Renumbered PDB file with chain_id_map does not match expected MD5 checksum"