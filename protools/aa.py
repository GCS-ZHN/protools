from Bio.Data import IUPACData
from typing import Dict, List
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


def validate_seq(s1: SeqLikeType, extra_symbols: str =''):
    """Validate sequence is only contain standard amino acid and allowed extra symbols."""
    s1 = ensure_seq_string(s1)
    aas = IUPACData.protein_letters + extra_symbols
    for pos, aa in enumerate(s1, 1):
        if aa not in aas:
            raise ValueError(f"Invalid amino acid: {aa} in sequence {s1} at {pos}")