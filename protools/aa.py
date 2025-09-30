from Bio.Data import IUPACData

RESIDUE2SIDECHAIN = {
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

RESIDUE2SIDECHAIN_3LETTER = {
    IUPACData.protein_letters_1to3[res]: atoms
    for res, atoms in RESIDUE2SIDECHAIN.items()
}