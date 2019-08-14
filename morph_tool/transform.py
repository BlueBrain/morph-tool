"""
Tools for morphology geometric transformations (translation, rotation, etc).
"""

import numpy as np


def _apply_recursively(func, obj, recenter=False):
    if hasattr(obj, 'soma'):
        obj.soma.points = func(obj.soma.points)
    for s in obj.iter():
        if recenter:
            s.points = s.points[0] + func(s.points - s.points[0])
        else:
            s.points = func(s.points)


def transform(obj, A):
    """
    Apply transformation matrix `A` to a given morphology object.

    Args:
        obj: Morphology / Section
        A: rotation matrix (4 x 4 NumPy array)
    """
    if A is None:
        return
    A = np.asarray(A)
    if A.shape != (4, 4):
        raise ValueError(
            "`A` should be 4 x 4 matrix (got instead: %s)" % str(A.shape)
        )
    A = A.transpose()
    func = lambda p: np.dot(np.column_stack((p, np.ones(len(p)))), A)[:, :3]
    _apply_recursively(func, obj)


def rotate(obj, A, recenter_section=False):
    """
    Apply rotation matrix `A` to a given morphology object.

    Args:
        obj: Morphology / Section
        A: rotation matrix (3 x 3 NumPy array)
        recenter_section (bool): recenter section before performing the rotation
            If False, the rotation is performed along the axis going through the origin
            If True, the rotation is performed along the axis going through the first point
                of each section
    """
    if A is None:
        return
    A = np.asarray(A)
    if A.shape != (3, 3):
        raise ValueError(
            "`A` should be 3 x 3 matrix (got instead: %s)" % str(A.shape)
        )
    A = A.transpose()
    func = lambda p: np.dot(p, A)
    _apply_recursively(func, obj, recenter=recenter_section)


def translate(obj, shift):
    """
    Apply translation to a given morphology object.

    Args:
        obj: Morphology / Section
        shift: shift vector ((x, y, z) NumPy array)
    """
    if shift is None:
        return
    shift = np.asarray(shift)
    if shift.shape != (3,):
        raise ValueError(
            "`shift` should be vector of shape (3,) (got instead: %s)" % str(shift.shape)
        )
    func = lambda p: p + shift
    _apply_recursively(func, obj)
