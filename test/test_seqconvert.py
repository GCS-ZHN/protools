from pathlib import Path
from protools import seqio
from protools import seqconvert
from .tools import md5_equal


MD5_DICT = dict()
with Path('data/fasta.md5').open() as f:
    for line in f:
        md5, filename = line.strip().split()
        filename = Path(filename)
        MD5_DICT[filename] = md5


def test_translate_from_rna():
    rna_fasta_file = Path("data/test_rna.fasta")
    rna_fasta = seqio.read_fasta(rna_fasta_file)
    protein_fasta_file = rna_fasta_file.with_suffix(".protein.standard.fasta")
    seqconvert.translate(
        rna_fasta,
        "Standard",
        nucleotide_type="rna").to_fasta(protein_fasta_file)
    assert md5_equal(protein_fasta_file, MD5_DICT[protein_fasta_file])

def test_translate_from_dna():
    dna_fasta_file = Path("data/test_dna.fasta")
    dna_fasta = seqio.read_fasta(dna_fasta_file)
    protein_fasta_file = dna_fasta_file.with_suffix(".protein.standard.fasta")
    seqconvert.translate(
        dna_fasta,
        "Standard",
        nucleotide_type="dna").to_fasta(protein_fasta_file)
    assert md5_equal(protein_fasta_file, MD5_DICT[protein_fasta_file])


def test_back_translate_to_rna():
    protein_fasta_file = Path("data/test_protein.fasta")
    protein_fasta = seqio.read_fasta(protein_fasta_file)
    rna_fasta_file = protein_fasta_file.with_suffix(".rna.standard.fasta")
    seqconvert.reverse_translate(
        protein_fasta,
        "Standard",
        nucleotide_type="rna").to_fasta(rna_fasta_file)
    assert md5_equal(rna_fasta_file, MD5_DICT[rna_fasta_file])


def test_back_translate_to_dna():
    protein_fasta_file = Path("data/test_protein.fasta")
    protein_fasta = seqio.read_fasta(protein_fasta_file)
    rna_fasta_file = protein_fasta_file.with_suffix(".dna.standard.fasta")
    seqconvert.reverse_translate(
        protein_fasta,
        "Standard",
        nucleotide_type="dna").to_fasta(rna_fasta_file)
    assert md5_equal(rna_fasta_file, MD5_DICT[rna_fasta_file])


def test_dna_optimizer():
    dna_fasta_file = Path("data/test_dna.fasta")
    dna_fasta = seqio.read_fasta(dna_fasta_file)
    optimized_dna_fasta_file = dna_fasta_file.with_suffix(".optimized.fasta")
    seqconvert.optimize_dna(
        dna_fasta,
        species="s_cerevisiae",
        avoid_patterns=["BsaI_site"]).to_fasta(optimized_dna_fasta_file)
