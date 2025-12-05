from Bio.Data import IUPACData
from typing import Callable, Dict, List, Literal
from protools.typedef import SeqLikeType
from protools.utils import ensure_seq_string


RESIDUE2SIDECHAIN: Dict[str, List[str]] = {
    'A': ['CB'],
    'R': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    'N': ['CB', 'CG', 'OD1', 'ND2'],
    'D': ['CB', 'CG', 'OD1', 'OD2'],
    'C': ['CB', 'SG'],
    'E': ['CB', 'CG', 'CD', 'OE1', 'OE2'],
    'Q': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
    'G': [],
    'H': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'I': ['CB', 'CG1', 'CG2', 'CD1'],
    'L': ['CB', 'CG', 'CD1', 'CD2'],
    'K': ['CB', 'CG', 'CD', 'CE', 'NZ'],
    'M': ['CB', 'CG', 'SD', 'CE'],
    'F': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'P': ['CB', 'CG', 'CD'],
    'S': ['CB', 'OG'],
    'T': ['CB', 'OG1', 'CG2'],
    'W': ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'Y': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
    'V': ['CB', 'CG1', 'CG2'],
}

RESIDUE2SIDECHAIN_3LETTER: Dict[str, List[str]] = {
    IUPACData.protein_letters_1to3[res]: atoms
    for res, atoms in RESIDUE2SIDECHAIN.items()
}


AA_PROPERTIES: Dict[str, List[str]] = {
    'hydrophobic': list('AILMFVPGW'),
    'positive': list('KRH'),
    'negative': list('DE'),
    'polar': list('STNQCY')
}

AA2PROPERTIES: Dict[str, str] = {
    aa: prop for prop, aas in AA_PROPERTIES.items() for aa in aas}


AAComparsionType = Literal['type', 'property'] | Callable[[str, str], bool]


def validate_seq(s1: SeqLikeType, extra_symbols: str =''):
    """Validate sequence is only contain standard amino acid and allowed extra symbols."""
    s1 = ensure_seq_string(s1)
    aas = IUPACData.protein_letters + extra_symbols
    for pos, aa in enumerate(s1, 1):
        if aa not in aas:
            raise ValueError(f"Invalid amino acid: {aa} in sequence {s1} at {pos}")


def aa_equal(a1: str, a2: str, comparsion: AAComparsionType = 'type') -> bool:
    """
    Compare two amino acids.

    Parameters
    ----------
        a1: str
            one-letter code of amino acid.
        a2: str
            one-letter code of amino acid.
        comparsion: string type of Callable[[str, str], bool]
            comparsion method for two amino acid. support bulitin methods
            'type' and 'property'. 'type' is the same amino acid, 'property' is the
            same amino acid property. If a callable is provided, it should take two
            amino acids and return a boolean indicating whether they are the same.
    
    Return
    ------
        bool
        Whether `a1` is equal to `a2` under provided comparsion criteria.

    Notes
    -----
        For non-standard or unknown aa token, this function assumes they 
        have different properties unless they are the same type tokens.
    
    """
    assert len(a1) == len(a2) == 1, "Only amino acid one letter code is allowed."
    a1 = a1.upper()
    a2 = a2.upper()
    if isinstance(comparsion, str):
        comparsion_funcs = {
            'type': lambda x, y: x==y,
            # assume unknown aa token have different property unless they are the same type
            'property': lambda x, y: (x==y) or (AA2PROPERTIES.get(x, 1) == AA2PROPERTIES.get(y, 2))
        }
        if comparsion in comparsion_funcs:
            comparsion = comparsion_funcs[comparsion]
        else:
            raise TypeError(f"Unknown comparsion type: {comparsion}")
    elif not isinstance(comparsion, Callable):
        raise TypeError(f"comparsion should be str or Callable type, not {type(comparsion)}")
    
    return comparsion(a1, a2)
