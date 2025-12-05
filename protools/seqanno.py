import warnings

from typing import Callable, Literal

from protools.utils import ensure_seq_string
from protools.typedef import SeqLikeType
from protools.aa import validate_seq, aa_equal, AAComparsionType
from functools import lru_cache


from antpack import SingleChainAnnotator, VJGeneTool
from Bio import Align


ChainType = Literal['H', 'L']
SchemeType = Literal['imgt', 'martin', 'kabat', 'aho']
SpecieType = Literal['human', 'mouse', 'alpaca', 'rabbit', 'unknown']
IdentityStrategyType = Literal[
    'min_aligned_length',
    'max_aligned_length',
    'avg_aligned_length',
    'min_sequence_length',
    'max_sequence_length',
    'first_alignment_length'
]

class SeqNumberingWarning(UserWarning):
    """
    Warning Defined for Seq Numbering.
    """
    pass


@lru_cache()
def get_annotator(chain: ChainType, scheme: SchemeType = 'imgt', is_tcr: bool = False) -> SingleChainAnnotator:
    """
    Get a SingleChainAnnotator for a given chain type.
    """
    _ANTIBODY_HC_IDS = ['H']
    _ANTIBODY_LC_IDS = ['K', 'L']
    _TCR_HC_IDS = ['B', 'D']
    _TCR_LC_IDS = ['A', 'G']
    _VALID_SCHEME = {
        'antibody': ['imgt', 'martin', 'kabat', 'aho'],
        'tcr': ['imgt']
    }
    if scheme == 'chothia':
        scheme = 'martin'
        warnings.warn(
            f"Chothia scheme is not supported, use Martin scheme (a modern version of Chothia) instead.",
            SeqNumberingWarning,
            stacklevel=2
            )

    seq_type = 'tcr' if is_tcr else 'antibody'
    if scheme not in _VALID_SCHEME[seq_type]:
        raise ValueError(
            f"Valid scheme for {seq_type} are {_VALID_SCHEME[seq_type]}, got {scheme}")

    if chain == 'H':
        annotator = SingleChainAnnotator(
            chains=_TCR_HC_IDS if is_tcr else _ANTIBODY_HC_IDS,
            scheme=scheme)
    else:
        annotator = SingleChainAnnotator(
            chains=_TCR_LC_IDS if is_tcr else _ANTIBODY_LC_IDS,
            scheme=scheme)
    
    return annotator


def numbering_seq(seq: SeqLikeType, chain: ChainType, is_tcr: bool = False, 
                  scheme: SchemeType = 'imgt') -> tuple[list[str], float, str, str]:
    """
    Numbering antibody or TCR sequences.

    Parameters
    ----------
    seq: str, Seq, or SeqRecord.
       The sequence to be numbered.

    chain: str
        'H' for heavy chain, 'L' for light chain
    
    is_tcr: bool
       True for TCR, False for antibody.

    scheme: str
        'imgt', 'kabat', 'aho' or 'martin' (modern version of chothia).

    Returns
    -------
    numbering position list
        A list of numbering positions
    sequence identity
        The sequence identity to the reference sequence.
    chain type
        The chain type, 'H', 'B', 'D' for heavy chain, 'K', 'L', 'A', 'G' for light chain.
    error message
        The error message if any.
    """
    seq = ensure_seq_string(seq)
    annotator = get_annotator(chain=chain, scheme=scheme, is_tcr=is_tcr)
    numbering, percent_identity, chain_type, err_message = annotator.analyze_seq(seq)
    if err_message:
        warnings.warn(err_message, stacklevel=2, category=SeqNumberingWarning)
    return numbering, percent_identity, chain_type, err_message


def anno_cdr(
        seq: SeqLikeType,
        chain: ChainType,
        is_tcr: bool = False,
        scheme: SchemeType = 'imgt',
        cdr_scheme: SchemeType| Literal['north'] = '') -> dict:
    """
    Annotate Antibody CDR regions based antpack package

    Parameters
    ----------
    seq : SeqLikeType
        Amino acid sequence of antibody

    chain : str
        Chain type, 'H' for heavy chain, 'K' for kappa light chain, 'L' for lambda light chain

    is_tcr: str
        Whether the sequence is TCR or not. Default is False. supported since antpack v0.3.8

    Returns
    ----------
    dict
        A dictionary containing CDR regions and their sequences
    """
    seq = ensure_seq_string(seq)
    numbering, percent_identity, anno_chain, _ = numbering_seq(
        seq, chain=chain, is_tcr=is_tcr, scheme=scheme)
    if len(numbering) != len(seq):
        raise ValueError(f"Numbering length {len(numbering)} != sequence length {len(seq)}")

    annotator = get_annotator(chain=chain, scheme=scheme, is_tcr=is_tcr)
    seq_range = {}
    label_map = lambda x: 'FR'+x[-1] if x.startswith('fmwk') else x.upper()
    cdr_labels = annotator.assign_cdr_labels(
        numbering=numbering, chain=anno_chain, scheme=cdr_scheme)
    for i , label in enumerate(cdr_labels):
        label = label_map(label)
        if label == '-':
            continue
        if label not in seq_range:
            seq_range[label] = [i, i+1]
        else:
            seq_range[label][1] = i + 1

    res = {}
    for label, seq_slice in seq_range.items():
        seq_slice = slice(*seq_slice)
        res[chain+label] = seq[seq_slice]
        res[chain+label+'_slice'] = seq_slice
    res[f'{chain}_percent_identity'] = percent_identity
    return res


def anno_vj_gene(seq: SeqLikeType,
                 chain_type: ChainType,
                 species: SpecieType,
                 is_tcr: bool = False) -> dict:
    """
    Annotate VJ gene for antibody or TCR sequence.
    
    Parameters
    ----------
    seq : SeqLikeType
        Amino acid sequence of antibody or TCR.
    
    chain_type : str
        Chain type, 'H' for heavy chain, 'K' for kappa light chain, 'L' for lambda light chain,
        'A' for alpha TCR chain, 'B' for beta TCR chain, 'D' for delta TCR chain, 'G' for gamma TCR chain.
    species : str
        Species of the sequence, can be 'human', 'mouse', 'alpaca', 'rabbit', or 'unknown'.
        if is_tcr is True, species must be 'human', 'mouse', or 'unknown'.
    is_tcr : bool
        Whether the sequence is TCR or not. Default is False.
    Returns
    -------
    dict
        A dictionary containing VJ gene information.
    """
    seq = ensure_seq_string(seq)
    vj_tool = VJGeneTool()
    alignment = numbering_seq(seq, chain_type, is_tcr=is_tcr)
    if is_tcr:
        assert species in ['human', 'mouse', 'unknown'], \
            "Species must be 'human', 'mouse' or 'unknown' for TCR"
    else:
        assert species in ['human', 'mouse', 'alpaca', 'rabbit', 'unknown'], \
            "Species must be 'human', 'mouse', 'alpaca', 'rabbit' or 'unknown' for antibody"
    v_gene, j_gene, v_pident, j_pident, species_matched = vj_tool.assign_vj_genes(
        alignment=alignment,
        sequence=seq,
        species=species,
        mode='identity'
    )
    if species != 'unknown':
        assert species_matched == species, \
            f"Species {species_matched} does not match the input species {species}"
    v_genes = v_gene.split('_')
    j_genes = j_gene.split('_')
    return {
        'v_genes': v_genes if len(v_genes) > 1  else v_genes[0],
        'j_genes': j_genes if len(j_genes) > 1 else j_genes[0],
        'v_pident': v_pident,
        'j_pident': j_pident,
        'species': species_matched,
    }


def calc_seq_identity(
        s1: SeqLikeType, s2: SeqLikeType, mode: Literal['global', 'local'] = 'global',
        strategy: IdentityStrategyType = 'min_aligned_length', max_alignments_size: int = 0) -> float:
    """
    Calculate sequence identity between two sequences.

    Parameters
    ----------
    s1 : SeqLikeType
        First sequence.
    s2 : SeqLikeType
        Second sequence.
    mode : str, optional
        Alignment mode, can be 'global' or 'local', by default 'global'.
    strategy : str, optional
        Strategy to calculate sequence identity, can be 'min_aligned_length', 'max_aligned_length',
        'avg_aligned_length', 'min_total_length', 'max_total_length', 'avg_total_length', 'min_sequence_length',
        'max_sequence_length', 'first_alignment_length',
        by default 'min_aligned_length'.
    max_alignments_size : int, optional
        Maximum number of alignments to consider, by default 0, which means all alignments.
        For large sequences, the number of alignments can be very large, which may cause
        memory overflow. In this case, you can set this parameter to a smaller value to limit
        the number of alignments to consider.
    Returns
    -------
    float
        Sequence identity between two sequences.
    """
    aligner = Align.PairwiseAligner(match_score = 1.0)
    if mode == 'global':
        aligner.mode = 'global'
    elif mode == 'local':
        aligner.mode = 'local'
    else:
        raise ValueError(f"Invalid mode: {mode}, must be 'global' or 'local'")
    
    alignments = aligner.align(s1, s2)
    try:
        alignments_size = len(alignments)
        if alignments_size == 0:
            return 0.0
        max_alignments_size = min(max_alignments_size, alignments_size) if max_alignments_size > 0 else alignments_size
        
    except OverflowError:
        if max_alignments_size == 0:
            max_alignments_size = 1000

    if strategy == 'min_aligned_length':
        length = min(alignment.length for i, alignment in enumerate(alignments) if i < max_alignments_size)
    elif strategy == 'max_aligned_length':
        length = max(alignment.length for i, alignment in enumerate(alignments) if i < max_alignments_size)
    elif strategy == 'avg_aligned_length':
        length = sum(alignment.length for i, alignment in enumerate(alignments) if i < max_alignments_size) / max_alignments_size
    elif strategy == 'min_sequence_length':
        length = min(len(s1), len(s2))
    elif strategy == 'max_sequence_length':
        length = max(len(s1), len(s2))
    elif strategy == 'first_alignment_length':
        length = alignments[0].length
    else:
        raise ValueError(f"Invalid strategy: {strategy}, must be 'min_aligned_length', 'max_aligned_length', 'avg_aligned_length', 'min_total_length', 'max_total_length', 'avg_total_length'")
    return alignments.score / length


def get_mutations(
        s1: SeqLikeType,
        s2: SeqLikeType,
        comparsion: AAComparsionType = 'type', chain_id: str = '') -> list[str]:
    """
    Generate mutations as a list of strings (e.g. ['A1W', 'C2D']).
    Position number (1-indexed) is the ungapped sequence position of `s1`.
    ['-2G', 'A2E', 'C4-'] means a 'G' insertion before 2rd aa in s1, a
    A-to-E replace at 2rd aa and a deletion at 4th aa.

    Parameters
    ----------
        s1: str, Seq or SeqRecord
            sequence of wt
        s2: str, Seq or SeqRecord
            sequnece of mutations, should be aligned with the same length 
            as wt.
        comparsion: string type of Callable[[str, str], bool]
            comparsion method for two amino acid. support bulitin methods
            'type' and 'property'. 'type' is the same amino acid, 'property' is the
            same amino acid property. If a callable is provided, it should take two
            amino acids and return a boolean indicating whether they are the same.
        chain_id: str
            chain id for the mutations, if provided, the mutations will be
            formatted as 'AH101Y'.

    Returns
    -------
       list of mutations in string format.
    """
    s1 = ensure_seq_string(s1)
    s2 = ensure_seq_string(s2)
    assert len(s1) == len(s2), "Sequences should be aligned with the same length."
    validate_seq(s1, extra_symbols='-')
    validate_seq(s2, extra_symbols='-')

    mutations = []

    alignment_pos_shift = 0
    for pos, (a1, a2) in enumerate(zip(s1, s2), 1):

        pos -= alignment_pos_shift
        if not aa_equal(a1, a2, comparsion):
            mutations.append(f'{a1}{chain_id}{pos}{a2}')

        if a1 == '-':
            alignment_pos_shift += 1
    return mutations
