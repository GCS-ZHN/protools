"""
A submodule struture operation.
"""

import numpy as np

from Bio.PDB.Entity import Entity
from scipy.spatial.transform import Rotation as R


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
