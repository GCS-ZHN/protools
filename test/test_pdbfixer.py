from protools import pdbfixer
from .tools import md5_equal


def test_renumber_residue():
    """
    Test the renumber_residue function.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        input_pdb = 'data/4ozi.pdb'
        exp_md5sum = '240b96faea97991fd0bcf747c0e7fe0a'
        output_pdb = f'{tmpdir}/4ozi_renumbered.pdb'
        pdbfixer.renumber_residue(
            pdb_file=input_pdb,
            out_file=output_pdb,
            model_idx=0,
            start=1
        )
        assert md5_equal(output_pdb, exp_md5sum), "Renumbered PDB file does not match expected MD5 checksum"


def test_renumber_residue_restart_by_chain():
    """
    Test the renumber_residue function with restart_by_chain option.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        input_pdb = 'data/4ozi.pdb'
        exp_md5sum = '612b2ecd5f1430a60c3987f110ea8378'  # Replace with actual expected MD5
        output_pdb = f'{tmpdir}/4ozi_renumbered_restart_by_chain.pdb'
        pdbfixer.renumber_residue(
            pdb_file=input_pdb,
            out_file=output_pdb,
            model_idx=0,
            start=1,
            restart_by_chain=True
        )
        assert md5_equal(output_pdb, exp_md5sum), "Renumbered PDB file with restart_by_chain does not match expected MD5 checksum"
