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
        entity: Entity, x_rotation: float = 0,
        y_rotation: float = 0, z_rotation: float = 0,
        degrees: bool = True,
        self_rotation: bool = False) -> np.ndarray:
    """
    Rotates the entity by the given angles inplace.

    Parameters
    ----------
    entity : Entity
        The entity to rotate.

    x_rotation : float
        The rotation angle around the x-axis.

    y_rotation : float
        The rotation angle around the y-axis.

    z_rotation : float
        The rotation angle around the z-axis.

    degrees : bool
        If True, the angles are interpreted as degrees, otherwise as radians.

    self_rotation : bool
        If True, the entity is rotated around its center of mass,
        otherwise around the origin.

    Returns
    -------
    np.ndarray
        The rotation matrix.
    """
    coord_point = np.array([0, 0, 0])
    if self_rotation:
        coord_point = np.mean([atom.get_coord() for atom in entity.get_atoms()], axis=0)
    rotation_matrix = R.from_euler(
        'xyz', [x_rotation, y_rotation, z_rotation], degrees=degrees).as_matrix()
    for atom in entity.get_atoms():
        coord = atom.get_coord() - coord_point
        new_coord = np.dot(rotation_matrix, coord) + coord_point
        atom.set_coord(new_coord)
    return rotation_matrix


def rand_rotate(
        entity: Entity, degrees: bool = True,
        self_rotation: bool = False) -> np.ndarray:
    """
    Rotates the entity by random angles.

    Parameters
    ----------
    entity : Entity
        The entity to rotate.

    degrees : bool
        If True, the angles are interpreted as degrees, otherwise as radians.

    self_rotation : bool
        If True, the entity is rotated around its center of mass, otherwise around the origin.

    Returns
    -------
    np.ndarray
        The rotation matrix.
    """
    x_rotation = np.random.uniform(0, 360)
    y_rotation = np.random.uniform(0, 360)
    z_rotation = np.random.uniform(0, 360)
    return rotate(entity, x_rotation, y_rotation, z_rotation, degrees, self_rotation)
