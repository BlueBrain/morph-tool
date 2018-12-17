"""
Tools for morphology geometric transformations (translation, rotation, etc).
"""

import numpy as np


def _apply_recursively(func, obj):
    if hasattr(obj, 'soma'):
        obj.soma.points = func(obj.soma.points)
    for s in obj.iter():
        s.points = func(s.points)


def transform(obj, A):
    """
    Apply transformation matrix `A` to a given morphology object.

    Args:
        obj: Morphology / Section
        A: rotation matrix (4 x 4 NumPy array)
    """
    A = np.asarray(A)
    if A.shape != (4, 4):
        raise ValueError(
            "`A` should be 4 x 4 matrix (got instead: %s)" % str(A.shape)
        )
    A = A.transpose()
    func = lambda p: np.dot(np.column_stack((p, np.ones(len(p)))), A)[:, :3]
    _apply_recursively(func, obj)


def rotate(obj, A):
    """
    Apply rotation matrix `A` to a given morphology object.

    Args:
        obj: Morphology / Section
        A: rotation matrix (3 x 3 NumPy array)
    """
    A = np.asarray(A)
    if A.shape != (3, 3):
        raise ValueError(
            "`A` should be 3 x 3 matrix (got instead: %s)" % str(A.shape)
        )
    A = A.transpose()
    func = lambda p: np.dot(p, A)
    _apply_recursively(func, obj)


def translate(obj, shift):
    """
    Apply translation to a given morphology object.

    Args:
        obj: Morphology / Section
        shift: shift vector ((x, y, z) NumPy array)
    """
    shift = np.asarray(shift)
    if shift.shape != (3,):
        raise ValueError(
            "`shift` should be vector of shape (3,) (got instead: %s)" % str(shift.shape)
        )
    func = lambda p: p + shift
    _apply_recursively(func, obj)
