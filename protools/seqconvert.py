"""
A module for converting between different sequence type.
"""

from .seqio import read_fasta, save_fasta
from Bio.Data import CodonTable
from pathlib import Path
from Bio.Seq import Seq


def _back_translate(seq: Seq, codon_table: CodonTable.CodonTable) -> Seq:
    """Back translate a protein sequence to a DNA/RNA sequence.

    Args:
        seq: Protein sequence.
        codon_table: Codon table to use for translation.

    Returns:
        DNA/RNA sequence.
    """
    result = ''
    for aa in seq:
        if aa == '*':
            result += '*'
        else:
            result += codon_table.back_table[aa]
    return Seq(result, codon_table.nucleotide_alphabet)



def translate(
        input_fasta: Path,
        output_fasta: Path,
        codon_table: str,
        nucleotide_type: str = 'dna'):
    """Translate a DNA/RNA fasta file to amino acid fasta file.

    Args:
        input_fasta: Path to input fasta file.
        output_fasta: Path to output fasta file.
        codon_table: Codon table to use for translation.
    """
    # Read input fasta file
    seqs = read_fasta(input_fasta)
    # Create codon table
    if nucleotide_type == 'rna':
        codon_table = CodonTable.unambiguous_rna_by_name[codon_table]
    elif nucleotide_type == 'dna':
        codon_table = CodonTable.unambiguous_dna_by_name[codon_table]
    else:
        raise ValueError(f'Invalid nucleotide type: {nucleotide_type}')
    # Translate sequences
    results = []
    for seq_record in seqs.values():
        seq_record.id = seq_record.id + '_translated'
        seq_record.seq = seq_record.seq.translate(table=codon_table, stop_symbol='')
        results.append(seq_record)
    save_fasta(results, output_fasta)


def back_translate(
        input_fasta: Path,
        output_fasta: Path,
        codon_table: str,
        nucleotide_type: str = 'dna'):
    """Reverse translate an amino acid fasta file to DNA/RNA fasta file.

    Args:
        input_fasta: Path to input fasta file.
        output_fasta: Path to output fasta file.
        codon_table: Codon table to use for translation.
    """
    # Read input fasta file
    seqs = read_fasta(input_fasta)
    # Create codon table
    if nucleotide_type == 'rna':
        codon_table = CodonTable.unambiguous_rna_by_name[codon_table]
    elif nucleotide_type == 'dna':
        codon_table = CodonTable.unambiguous_dna_by_name[codon_table]
    else:
        raise ValueError(f'Invalid nucleotide type: {nucleotide_type}')
    # Translate sequences
    results = []
    
    for seq_record in seqs.values():
        seq_record.id = seq_record.id + '_reverse_translated'
        seq_record.seq = _back_translate(seq_record.seq, codon_table)
        results.append(seq_record)
    save_fasta(results, output_fasta)
