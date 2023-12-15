from pathlib import Path
from protools import seqconvert
from .tools import md5_equal


def test_translate_from_rna():
    rna_fasta = Path("data/test_rna.fasta")
    protein_fasta = rna_fasta.with_suffix(".protein.standard.fasta")
    seqconvert.translate(
        rna_fasta,
        protein_fasta,
        "Standard",
        nucleotide_type="rna")
    expected_md5 = "84ab7b2a20cad2e319fca14bae504f8a"
    assert md5_equal(protein_fasta, expected_md5)
    protein_fasta = rna_fasta.with_suffix(".protein.yeast.fasta")
    seqconvert.translate(
        rna_fasta,
        protein_fasta,
        "Yeast Mitochondrial",
        nucleotide_type="rna")
    expected_md5 = "2abf78a63dd9fdb3ab47d020f60b08c4"
    assert md5_equal(protein_fasta, expected_md5)


def test_translate_from_dna():
    rna_fasta = Path("data/test_dna.fasta")
    protein_fasta = rna_fasta.with_suffix(".protein.standard.fasta")
    seqconvert.translate(
        rna_fasta,
        protein_fasta,
        "Standard",
        nucleotide_type="dna")
    expected_md5 = "fa6f3ea2e2f4201c2513669e17a6650d"
    assert md5_equal(protein_fasta, expected_md5)
    protein_fasta = rna_fasta.with_suffix(".protein.yeast.fasta")
    seqconvert.translate(
        rna_fasta,
        protein_fasta,
        "Yeast Mitochondrial",
        nucleotide_type="dna")
    expected_md5 = "1734be71daa8ec519a409d8deebeb718"
    assert md5_equal(protein_fasta, expected_md5)


def test_back_translate_to_rna():
    protein_fasta = Path("data/test_protein.fasta")
    rna_fasta = protein_fasta.with_suffix(".rna.standard.fasta")
    seqconvert.back_translate(
        protein_fasta,
        rna_fasta,
        "Standard",
        nucleotide_type="rna")
    expected_md5 = "e73aa4dde269996948e4221a6d86536f"
    assert md5_equal(rna_fasta, expected_md5)
    rna_fasta = protein_fasta.with_suffix(".rna.yeast.fasta")
    seqconvert.back_translate(
        protein_fasta,
        rna_fasta,
        "Yeast Mitochondrial",
        nucleotide_type="rna")
    expected_md5 = "93aa63ba0c15f1ba7f2d0ef241471828"
    assert md5_equal(rna_fasta, expected_md5)


def test_back_translate_to_dna():
    protein_fasta = Path("data/test_protein.fasta")
    rna_fasta = protein_fasta.with_suffix(".dna.standard.fasta")
    seqconvert.back_translate(
        protein_fasta,
        rna_fasta,
        "Standard",
        nucleotide_type="dna")
    expected_md5 = "2976f4323e08c2d8038d2a1564573ec7"
    assert md5_equal(rna_fasta, expected_md5)
    rna_fasta = protein_fasta.with_suffix(".dna.yeast.fasta")
    seqconvert.back_translate(
        protein_fasta,
        rna_fasta,
        "Yeast Mitochondrial",
        nucleotide_type="dna")
    expected_md5 = "ee778dc154934e476da2a8db640a7be0"
    assert md5_equal(rna_fasta, expected_md5)
