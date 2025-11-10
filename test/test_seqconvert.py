import pytest

from protools import seqio
from protools import seqconvert
from .tools import md5_equal


@pytest.mark.parametrize(
        'seq_file, nucleotide_type, exp_md5sum',
        [
            ('data/test_rna.fasta', 'rna', 'c3d53f7293115ea301378ee28afbdd2b'),
            ('data/test_dna.fasta', 'dna', '093e3f40428c88f7629f1dcfe3bd1b23')
         ]
)
def test_translate(tmp_path, seq_file, nucleotide_type, exp_md5sum):
    rna_fasta = seqio.read_fasta(seq_file)
    protein_fasta_file = tmp_path / "translated_protein.fasta"
    seqconvert.translate(
        rna_fasta,
        "Standard",
        nucleotide_type=nucleotide_type).to_fasta(protein_fasta_file)
    assert md5_equal(protein_fasta_file, exp_md5sum)


@pytest.mark.parametrize(
        'seq_file, nucleotide_type, exp_md5sum',
        [
            ('data/test_protein.fasta', 'rna', '5573fc2e28931803b6b7bacaa576a004'),
            ('data/test_protein.fasta', 'dna', 'f905beb77901a82fb4628bc6015a0374')
         ]
)
def test_reverse_translate(tmp_path, seq_file, nucleotide_type, exp_md5sum):
    protein_fasta = seqio.read_fasta(seq_file)
    rna_fasta_file = tmp_path / "reverse_translated_nucleotide.fasta"
    seqconvert.reverse_translate(
        protein_fasta,
        "Standard",
        nucleotide_type=nucleotide_type).to_fasta(rna_fasta_file)
    assert md5_equal(rna_fasta_file, exp_md5sum)


@pytest.mark.parametrize(
        'seq_file, species, avoid_patterns',
        [
            ('data/test_dna.fasta', 's_cerevisiae', ['BsaI_site'])
        ]
)
def test_dna_optimizer(tmp_path, seq_file, species, avoid_patterns):
    dna_fasta = seqio.read_fasta(seq_file)
    optimized_dna_fasta_file = tmp_path / "dna_optimized.fasta"
    seqconvert.optimize_dna(
        dna_fasta,
        species=species,
        avoid_patterns=avoid_patterns).to_fasta(optimized_dna_fasta_file)
