"""
Tools for morphology geometric transformations (translation, rotation, etc).
"""
import logging
from enum import Enum
import numpy as np

from scipy.spatial.transform import Rotation

from morphio import SectionType, IterType
from neurom.morphmath import angle_between_vectors

from morph_tool.spatial import point_to_section_segment
from morph_tool.apical_point import apical_point_section_segment

L = logging.getLogger(__name__)


def _apply_recursively(func, obj, origin=(0, 0, 0)):
    origin = np.array(origin)

    if hasattr(obj, 'soma'):
        obj.soma.points = origin + func(obj.soma.points - origin)
    for s in obj.iter():
        s.points = origin + func(s.points - origin)


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


def rotate(obj, A, origin=(0, 0, 0)):
    """
    Apply rotation matrix `A` to a given morphology object.

    Args:
        obj: Morphology / Section
        A: rotation matrix (3 x 3 NumPy array)
        origin (3D point): the origin of the rotation
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
    _apply_recursively(func, obj, origin)


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


def align(section, direction):
    '''Rotate a section (and all its descendents) so that its initial segment is oriented along
    "direction"'''
    section_dir = section.points[1] - section.points[0]
    alpha = angle_between_vectors(section_dir, direction)
    if alpha < 1e-8:
        return

    if abs(alpha - np.pi) < 1e-8:
        axis = np.cross(section_dir, [1, 0, 0])

        # Case where X axis and section_dir are colinear
        if np.linalg.norm(axis) < 1e-8:
            axis = np.cross(section_dir, [0, 1, 0])
    else:
        axis = np.cross(section_dir, direction)
    axis /= np.linalg.norm(axis)
    matrix = Rotation.from_rotvec(alpha * axis).as_matrix()

    rotate(section, matrix, origin=section.points[0])


class AlignMethod(Enum):
    """Contains possible align methods for align_morphology"""

    WHOLE = 'whole'
    TRUNK = 'trunk'
    FIRST_SECTION = 'first_section'
    FIRST_SEGMENT = 'first_segment'

    @classmethod
    def values(cls):
        """Get all possible values."""
        return list(map(lambda c: c.value, cls))


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2

    Picked from: https://stackoverflow.com/a/59204638/3868743
    Args:
        vec1: A 3d "source" vector
        vec2: A 3d "destination" vector

    Returns:
        A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    vec1, vec2 = vec1 / np.linalg.norm(vec1), vec2 / np.linalg.norm(vec2)

    v_cross = np.cross(vec1, vec2)
    v_cross_norm = np.linalg.norm(v_cross)
    if v_cross_norm == 0:
        return np.eye(3)

    kmat = np.array([[0.0, -v_cross[2], v_cross[1]],
                     [v_cross[2], 0.0, -v_cross[0]],
                     [-v_cross[1], v_cross[0], 0.0]])

    return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - np.dot(vec1, vec2)) / (v_cross_norm ** 2))


# pylint: disable=inconsistent-return-statements
def _get_points(morph, method, neurite_type, target_point):
    """Extract relevant points of dendrite to align the morphology, see align_morphology."""
    _to_type = {'apical': SectionType.apical_dendrite, 'axon': SectionType.axon}

    for root_section in morph.root_sections:
        if root_section.type == _to_type[neurite_type]:
            if method == AlignMethod.TRUNK.value:
                if target_point is not None:
                    target_secid = point_to_section_segment(morph, target_point)[0] - 1
                    if target_secid is None:
                        return None
                elif neurite_type == 'apical':
                    target_secid = apical_point_section_segment(morph)[0]
                else:
                    raise Exception(f"We don't know how to get target point for {neurite_type}.")

                return np.vstack(
                    [section.points
                     for section in morph.sections[target_secid].iter(IterType.upstream)]
                )
            if method == AlignMethod.FIRST_SECTION.value:
                return root_section.points

            if method == AlignMethod.FIRST_SEGMENT.value:
                return root_section.points[:2]

            return np.vstack([section.points for section in root_section.iter()])


def _get_principal_direction(points):
    '''Return the principal direction of a point cloud
    It is the eigen vector of the covariance matrix with the highest eigen value.

    Taken from neuror.unravel.'''
    X = np.copy(np.asarray(points))
    X -= np.mean(X, axis=0)
    C = np.dot(X.T, X)
    w, v = np.linalg.eig(C)
    return v[:, w.argmax()]


def align_morphology(
    morph, direction=None, method='whole', neurite_type='apical', target_point=None
):
    """In-place alignment of a morphology towards a 'direction'.

    The base algorithm is based on eigenvalue decomposition of the correlation matrix obtained
    from points in specified neurites, giving the principal axis of the neurite.
    Currently, five algorithms are implemented, differing in the choice of points:

    1) with method='whole': All the points in the apical dendrite are used.
    2) with method='trunk': Points in section up to the target points are used.
    3) with method='first_section': Only the points in the first section are used.
    4) with method='first_segment': Only the points in the first segment are used.
    5) with method an ndarray or list, we will use it as the direction directly

    If no neurite is present, no rotation is applied, and the identity rotation is returned.
    If two neurites of same types are present, the first accessed by Morphio will be used.

    Args:
        morph (morphio.Morphology): morphology to align
        direction (ndarray): 3-vector for final direction, if None, [0, 1, 0] will be used
        method (str|ndarray): method for alignment.
        neurite_type (str): neurite to consider, can only be apical or axon
        target_point (ndarray): position of target point for method='trunk',
            if None and neurite_type='apical', it will be estimated

    Returns:
        3x3 array with applied rotation matrix
    """
    if isinstance(method, str) and method not in AlignMethod.values():
        raise NotImplementedError(f"Method {method} is not implementd")

    if direction is None:
        direction = [0.0, 1.0, 0.0]
    else:
        direction /= np.linalg.norm(direction)

    if isinstance(method, (np.ndarray, list)):
        points = np.array([[0.0, 0.0, 0.0], method])
    else:
        points = _get_points(morph, method, neurite_type, target_point)

    if points is None:
        L.info('We did not find an apical point to align the morphology')
        return np.eye(3)

    principal_direction = _get_principal_direction(points)
    principal_direction *= np.sign(points.dot(principal_direction).sum())

    rotation_matrix = rotation_matrix_from_vectors(principal_direction, direction)
    rotate(morph, rotation_matrix)

    return rotation_matrix
