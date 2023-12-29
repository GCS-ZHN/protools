"""
A module for converting between different sequence type.
"""

from .seqio import Fasta, read_fasta, save_fasta
from Bio.Data import CodonTable
from Bio.Seq import Seq
from .typedef import FilePathType
from .utils import ensure_path
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceGCContent, EnforceTranslation, CodonOptimize


def _reverse_translate(seq: Seq, codon_table: CodonTable.CodonTable) -> Seq:
    """Back translate a protein sequence to a DNA/RNA sequence.

    Parameters
    ----------
    seq : Seq
        Protein sequence to back translate.
    codon_table : CodonTable.CodonTable
        Codon table to use for back translation.

    Returns
    ----------
    Seq
        Back translated DNA/RNA sequence.
    """
    result = ''
    for aa in seq:
        if aa == '*':
            result += '*'
        else:
            result += codon_table.back_table[aa]
    return Seq(result.lower())


def translate(
        input_fasta: Fasta,
        codon_table: str = 'Standard',
        nucleotide_type: str = 'dna', 
        inplace: bool = False) -> Fasta:
    """Translate a DNA/RNA fasta file to amino acid fasta file.

    Parameters
    ----------
    input_fasta : Fasta
        Input fasta object.
    codon_table : str
        Codon table to use for translation.
    nucleotide_type : str
        Type of nucleotide to translate.
    inplace : bool
        Whether to modify the input fasta object in place.

    Returns
    ----------
    Fasta
        Translated fasta object.
    """
    # Create codon table
    if nucleotide_type == 'rna':
        codon_table = CodonTable.unambiguous_rna_by_name[codon_table]
    elif nucleotide_type == 'dna':
        codon_table = CodonTable.unambiguous_dna_by_name[codon_table]
    else:
        raise ValueError(f'Invalid nucleotide type: {nucleotide_type}')
    # Translate sequences
    results = input_fasta if inplace else Fasta()
    for seq_record in input_fasta.values():
        if not inplace:
            seq_record = SeqRecord(
                seq=seq_record.seq,
                id=seq_record.id,
                description='',
                name='')
        seq_record.id = seq_record.id + '_translated'
        seq_record.description = ''
        seq_record.seq = seq_record.seq.translate(table=codon_table, stop_symbol='')
        if not inplace:
            results[seq_record.id] = seq_record

    return results


def reverse_translate(
        input_fasta: Fasta,
        codon_table: str = 'Standard',
        nucleotide_type: str = 'dna',
        inplace: bool = False) -> Fasta:
    """Reverse translate an amino acid fasta file to DNA/RNA fasta file.

    Parameters
    ----------
    input_fasta : Fasta
        Input fasta object.
    codon_table : str
        Codon table to use for translation.
    nucleotide_type : str
        Type of nucleotide to translate.
    inplace : bool
        Whether to modify the input fasta object in place.

    Returns
    ----------
    Fasta
        Translated fasta object.
    """
    # Create codon table
    if nucleotide_type == 'rna':
        codon_table = CodonTable.unambiguous_rna_by_name[codon_table]
    elif nucleotide_type == 'dna':
        codon_table = CodonTable.unambiguous_dna_by_name[codon_table]
    else:
        raise ValueError(f'Invalid nucleotide type: {nucleotide_type}')
    # Translate sequences
    results = input_fasta if inplace else Fasta()
    for seq_record in input_fasta.values():
        if not inplace:
            seq_record = SeqRecord(
                seq=seq_record.seq,
                id=seq_record.id,
                description='',
                name='')
        seq_record.id = seq_record.id + '_reverse_translated'
        seq_record.description = ''
        seq_record.seq = _reverse_translate(seq_record.seq, codon_table)
        if not inplace:
            results[seq_record.id] = seq_record
    return results


def optimize_dna(
        input_fasta: Fasta,
        species: str,
        codon_table: str = 'Standard',
        avoid_patterns: list = None,
        inplace: bool = False) -> Fasta:
    """Optimize a DNA fasta file to maximize expression.

    Parameters
    ----------
    input_fasta : Fasta
        Input fasta object.
    species : str
        Species to optimize for.
    codon_table : str
        Codon table to use for translation.
    avoid_patterns : list
        List of patterns to avoid.
    inplace : bool
        Whether to modify the input fasta object in place.

    Returns
    ----------
    Fasta
        Optimized fasta object.
    """
    # Read input fasta file
    results = input_fasta if inplace else Fasta()
    constraints = [
        EnforceGCContent(mini=0.3, maxi=0.7, window=50),
        EnforceTranslation(genetic_table=codon_table)
    ]
    if avoid_patterns is not None:
        constraints += [AvoidPattern(pattern=pattern) for pattern in avoid_patterns]
    for seq_record in input_fasta.values():
        if not inplace:
            seq_record = SeqRecord(
                seq=seq_record.seq,
                id=seq_record.id,
                description='',
                name='')
        seq_record.id = seq_record.id + '_optimized'
        seq_record.description = ''
        problem = DnaOptimizationProblem(
            sequence=str(seq_record.seq),
            constraints=constraints,
            objectives=[CodonOptimize(species=species)]
        )
        problem.resolve_constraints()
        new_seq = Seq(problem.sequence.lower())
        assert new_seq.translate(table=codon_table) == seq_record.seq.translate(table=codon_table), \
            f'Optimization failed for {seq_record.id}'
        seq_record.seq = new_seq
        if not inplace:
            results[seq_record.id] = seq_record
    return results


class ComplexFasta(OrderedDict):
    """
    A class for handling fasta files with multiple sequences per record.
    """
        
    @staticmethod
    def read_colabfold(source: FilePathType) -> 'ComplexFasta':
        """Read a fasta file in colabfold format.

        Parameters
        ----------
        source : FilePathType
            Path to fasta file.

        Returns
        ----------
        ComplexFasta
            Fasta object.
        """
        source = ensure_path(source)
        results = ComplexFasta()
        raw_fasta = read_fasta(source)
        for seq_id, seq_record in raw_fasta.items():
            complex_seqs = str(seq_record.seq).split(':')
            results[seq_id] = Fasta((f'{seq_id}_{i}', seq) for i, seq in enumerate(complex_seqs))
        return results

    @staticmethod
    def read_rosettafold(source: FilePathType) -> 'ComplexFasta':
        """Read a fasta file in rosettafold format.

        Parameters
        ----------
        source : FilePathType
            Path to fasta directory.

        Returns
        ----------
        ComplexFasta
            Fasta object.
        """
        source = ensure_path(source)
        if not source.is_dir():
            raise ValueError(f'{source} is not a directory.')
        results = ComplexFasta()
        for fasta_file in source.glob('*.fasta'):
            seq_id = fasta_file.stem
            results[seq_id] = read_fasta(fasta_file)
        return results
        
    def _write_by_sep(self, target: FilePathType, sep: str) -> None:
        records = (SeqRecord(
            seq=Seq(sep.join(str(seq_record.seq) for seq_record in fasta.values())),
            id=seq_id,
            description='',
            name=''
        ) for seq_id, fasta in self.items())
        save_fasta(records, target)

    def write_colabfold(self, target: FilePathType) -> None:
        """Write complex in colabfold format.

        Parameters
        ----------
        target : FilePathType
            Path to fasta file.
        """
        self._write_by_sep(target, ':')

    def write_af2rank(self, target: FilePathType) -> None:
        """Write complex in af2rank format.

        Parameters
        ----------
        target : FilePathType
            Path to fasta file.
        """
        self._write_by_sep(target, '')

    def write_rosettafold(self, target: FilePathType) -> None:
        """Write complex in rosettafold format.

        Parameters
        ----------
        target : FilePathType
            Path to fasta directory.
        """
        target = ensure_path(target)
        target.mkdir(parents=True)
        for seq_id, fasta in self.items():
            fasta.to_fasta(target / f'{seq_id}.fasta')
    
    def __setitem__(self, __key: str, __value: Fasta) -> None:
        if not isinstance(__value, Fasta):
            raise ValueError(f'Value must be a Fasta object, got {type(__value)}')
        if __key in self:
            raise ValueError(f'Key {__key} already exists')
        return super().__setitem__(__key, __value)
    
    def __getitem__(self, __key: str) -> Fasta:
        return super().__getitem__(__key)