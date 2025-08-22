"""
A submodule struture operation.
"""

import numpy as np
import re

from Bio.PDB.Entity import Entity
from Bio.PDB.Chain import Chain
from scipy.spatial.transform import Rotation as R

from protools import pdbio, utils

from typing import List


__all__ = ['translate', 'rotate', 'rand_rotate']


def translate(
        entity: Entity, x_translation: float = 0,
        y_translation: float = 0, z_translation: float = 0) -> np.ndarray:
    """
    Translates the entity by the given vector inplace.

    Parameters
    ----------
    entity : Entity
        The entity to translate.

    x_translation : float
        The translation along the x-axis.

    y_translation : float
        The translation along the y-axis.

    z_translation : float
        The translation along the z-axis.

    Returns
    -------
    np.ndarray
        The translation vector.
    """
    translation_vector = np.array([x_translation, y_translation, z_translation])
    for atom in entity.get_atoms():
        coord = atom.get_coord()
        new_coord = coord + translation_vector
        atom.set_coord(new_coord)
    return translation_vector


def rotate(
        entity: Entity, 
        rotation_angles: np.ndarray,
        degrees: bool = True,
        self_rotation: bool = False) -> Entity:
    """
    Rotates the entity by the given angles inplace.

    Parameters
    ----------
    entity : Entity
        The entity to rotate.

    rotation_angles : np.ndarray

    degrees : bool
        If True, the angles are interpreted as degrees, otherwise as radians.

    self_rotation : bool
        If True, the entity is rotated around its center of mass,
        otherwise around the origin.

    Returns
    -------
    Entity
        The rotated entity (inplace).
    """
    rotation_angles = np.array(rotation_angles)
    if len(rotation_angles) != 3 or len(rotation_angles.shape) != 1:
        raise ValueError('The rotation angles must be a 1D array of length 3.')
    coord_point = np.array([0, 0, 0])
    if self_rotation:
        coord_point = np.mean(
            [atom.get_coord() for atom in entity.get_atoms()], axis=0)
    rotation_matrix = R.from_euler(
        'xyz', rotation_angles, degrees=degrees).as_matrix()
    for atom in entity.get_atoms():
        coord = atom.get_coord() - coord_point
        new_coord = np.dot(rotation_matrix, coord) + coord_point
        atom.set_coord(new_coord)
    return entity

def rand_rotate(
        entity: Entity, degrees: bool = True,
        self_rotation: bool = False,
        seed: int = None) -> np.ndarray:
    """
    Rotates the entity by random angles.

    Parameters
    ----------
    entity : Entity
        The entity to rotate.

    degrees : bool
        If True, the angles are interpreted as degrees, otherwise as radians.

    self_rotation : bool
        If True, the entity is rotated around its center of mass,
        otherwise around the origin.

    seed : int
        The seed for the random number generator.

    Returns
    -------
    np.ndarray
        The rotation angles.
    """
    rand_gen = np.random.default_rng(seed)
    if degrees:
        rotation_angles = rand_gen.uniform(0, 360, 3)
    else:
        rotation_angles = rand_gen.uniform(0, 2 * np.pi, 3)
    rotate(entity, rotation_angles, degrees, self_rotation)
    return rotation_angles


def chain_split(
        chain: Chain,
        linker_pattern: str) -> List[Chain]:
    """
    Split chain by the given regular expression linker pattern.

    Parameters
    ----------
    chain : Chain
        The chain to split.

    linker_pattern : str
        The regular expression linker pattern. such as 
        'GGGGS', 'G{3,4}S'
    """

    aa_sequence = str(pdbio.get_aa_sequence(chain))
    sequence_size = len(aa_sequence)
    aa_residues = list(pdbio.get_aa_residues(chain).values())
    parent: Entity = chain.get_parent()
    if parent is not None:
        parent.detach_child(chain.get_id())
    assert sequence_size == len(aa_residues)

    pattern = re.compile(linker_pattern)
    # found all linker slices
    matches = pattern.finditer(aa_sequence)
    linker_slices = [slice(*m.span()) for m in matches]
    linker_intervals = utils.Intervals.from_slices(linker_slices)
    new_chain_intervals = linker_intervals.invert(high=sequence_size)
    new_chains = []
    idx = 0
    for interval in new_chain_intervals:
        while True:
            if idx > 25:
                raise ValueError('Too many chains splited.')
            chain_id = chr(ord('A') + idx)
            if parent is not None and parent.has_id(chain_id):
                idx += 1
                continue
            break
        new_chain = Chain(chain_id)
        for residue in aa_residues[interval]:
            new_chain.add(residue)
        new_chains.append(new_chain)
        idx += 1
        if parent is not None:
            parent.add(new_chain)
    return new_chains
