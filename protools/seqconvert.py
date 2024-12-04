"""
A module for converting between different sequence type.
"""
import numpy as np
import pandas as pd

from .seqio import Fasta, read_fasta, save_fasta
from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq
from .typedef import FilePathType
from .utils import ensure_path
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceGCContent, EnforceTranslation, CodonOptimize
from protools.seqanno import numbering_seq
from typing import Iterable


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


class AntibodyPSSM(object):
    """
    A class for building antibody PSSM.

    Parameters
    ----------
    pssm : np.ndarray
        The PSSM matrix.
    """

    def __init__(self, pssm: np.ndarray = None):
        self._pssm = pssm

    @classmethod
    def build_from_fasta(cls, fasta: Fasta, chain: str = None) -> 'AntibodyPSSM':
        """Build PSSM from a fasta object of antibody sequences.

        Parameters
        ----------
        fasta : Fasta
            Input fasta object containing antibody sequences.
        chain : str
            Chain type, 'H' for heavy chain, 'L' for lambda/kappa light chain.
            if sequence already aligned, set chain to None.

        Returns
        ----------
        AntibodyPSSM
            An instance of AntibodyPSSM with the built PSSM.
        """
        cls.build_from_seqs(map(str, fasta.values()), chain)

    @classmethod
    def build_from_seqs(cls, seqs: Iterable[str], chain: str = None) -> 'AntibodyPSSM':
        """Build PSSM from a list of antibody sequences.

        Parameters
        ----------
        seqs : Iterable[str]
            Input sequences.
        chain : str
            Chain type, 'H' for heavy chain, 'L' for lambda/kappa light chain.
            if sequence already aligned, set chain to None.
        """
        # Initialize a dictionary to hold frequency counts for each position
        position_counts = {}

        # Iterate over each sequence in the fasta
        aligned_length = None
        for seq in seqs:
            # Get numbering using the numbering_seq function
            if chain is not None:
                numbering, _ = numbering_seq(seq, chain)
            else:
                seq_length = len(seq)
                if aligned_length is None:
                    aligned_length = seq_length
                if seq_length != aligned_length:
                    raise ValueError(
                        "Aligned sequence length not matched."
                    )
                numbering = [i if s in IUPACData.protein_letters else '-' for i, s in enumerate(seq)]

            # Count amino acid frequencies at each position
            for seq_idx, num in enumerate(numbering):
                try:
                    # Skip non-Fv region amino acids and insertions
                    num = int(num)
                except ValueError:
                    continue
                if num not in position_counts:
                    position_counts[num] = {aa: 0 for aa in IUPACData.protein_letters}
                position_counts[num][seq[seq_idx]] += 1

        # Convert counts to frequencies
        pssm = []
        for position in sorted(position_counts.keys()):
            total = sum(position_counts[position].values())
            frequencies = [position_counts[position][aa] / total for aa in IUPACData.protein_letters]
            pssm.append(frequencies)

        # Convert the list of frequencies to a numpy array
        pssm_array = np.array(pssm)
        pssm_array = np.log(pssm_array * 20 + 1e-12)
        # Return an instance of AntibodyPSSM
        return cls(pssm=pssm_array)

    @classmethod
    def read(cls, source: FilePathType) -> 'AntibodyPSSM':
        with open(source, 'rb') as f:
            pssm = np.load(f)
        return cls(pssm)

    def write(self, target: FilePathType):
        with open(target, 'wb') as f:
            np.save(f, self.pssm)

    @property
    def pssm(self) -> np.ndarray:
        """The PSSM matrix."""
        if self._pssm is None:
            raise ValueError('PSSM not built')
        return self._pssm

    @property
    def logits(self) -> np.ndarray:
        """The logit matrix."""
        return self.pssm - np.log(20)
    
    @property
    def probs(self) -> np.ndarray:
        """The protein probability matrix (PPM)."""
        e_x = np.exp(self.logits)
        return e_x / e_x.sum(axis=1, keepdims=True)
    
    def __repr__(self) -> str:
        if self.pssm is None:
            return 'AntibodyPSSM(None)'
        
        pssm_df = pd.DataFrame(self.pssm, columns=list(IUPACData.protein_letters))
        return repr(pssm_df)
    
    def sample(self, seed: int = None) -> str:
        """Sample a sequence from the PSSM.

        Parameters
        ----------
        seed : int
            Random seed.

        Returns
        ----------
        str
            Sampled sequence.
        """
        if seed is not None:
            np.random.seed(seed)
        seq = ''
        for row in self.probs:
            seq += np.random.choice(list(IUPACData.protein_letters), p=row)
        return seq


def clip_a3m(a3m: Fasta, start: int = 0, end: int = None) -> Fasta:
    """
    Clip single a3m file


    Parameters
    ---------------
    - start: int, default is 0
            the inclusive beginning of clipped sequence.
    
    - end: int, default is None
            the exclusive ending of clipped sequence.
    """
    items = iter(a3m.items())
    query_id, query_seq = next(items)
    query_seq = str(query_seq.seq)
    if end is None:
        end = len(query_seq)
    if end > len(query_seq) or end <= start or start < 0:
        raise ValueError('Invalid position for clip.')
    
    new_fasta = Fasta(query_id=query_seq[start:end])

    for rid, seq in items:
        seq = str(seq.seq)
        s, e = 0, 0
        aligned_idx = -1
        for idx, aa in enumerate(seq):
            if aa.isupper() or aa == '-':
                aligned_idx += 1
                if aligned_idx == start:
                    s = idx
                if aligned_idx == end - 1:
                    e = idx + 1
                    break
            elif not aa.islower():
                raise ValueError(f'Invalid symol for seq {rid} at position {idx+1}')
        else:
            raise RuntimeError(f'Clip failed for {rid}')
        new_fasta[rid] = seq[s:e]

    return new_fasta
