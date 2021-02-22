"""
Tools for morphology geometric transformations (translation, rotation, etc).
"""
import logging
import numpy as np

from scipy.spatial.transform import Rotation
from scipy.linalg import svd

from morphio import SectionType, IterType
from neurom.morphmath import angle_between_vectors

from morph_tool.apical_point import apical_point_section_segment
from morph_tool.spatial import point_to_section_segment

logger = logging.getLogger('morph_tool')


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


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2

    Picked from: https://stackoverflow.com/a/59204638/3868743

    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s == 0:
        return np.eye(3)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


# pylint: disable=too-many-branches
def align_pyr_morph(morph, direction=None, method='whole_apical', apical_point=None):
    """In-place alignment of a pyramidal morphology such that apical tree is along 'direction'.

    The base algorithm is based on SVD decomposition of the correlation matrix obtained from
    points in apical neurites, giving the principal axis of apical dendrite.
    Currently, three algorithms are implemented, differing in the choice of points:

    1) with method='whole_apical': All the points in the apical dendrite are used.
    2) with method='first_segment': Only the 2 first points in the first section are used.
    3) with method='apical_trunk': Points in section up to the apical points are used.

    If no apical is present, no rotation is applied, and the identity rotation is returned.
    If two apicals are present, the apical accessed in second by Morphio will be used.

    Args:
        morph (morphio.Morphology): morphology to align
        direction (ndarray): 3-vector for final direction, if None, [0, 1, 0] will be used
        method (str): method for alignment.
        apical_point (ndarray): position of apical point for method='apical_trunk',
            if None, it will be estimated

    Returns:
        3x3 array with applied rotation matrix
    """
    if method not in ['whole_apical', 'apical_trunk', 'first_segment']:
        raise NotImplementedError(f'Method {method} is not implemented')

    if direction is None:
        direction = [0.0, 1.0, 0.0]
    apical_points = None
    apical_root_sections = [
        sec for sec in morph.root_sections if sec.type == SectionType.apical_dendrite
    ]

    if len(apical_root_sections) > 1:
        logger.warning(
            'The morphology has %s apical dendrites, only the first one will be used.',
            len(apical_root_sections)
        )
    elif len(apical_root_sections) == 0:
        logger.info('No apical dendrite was found in the morphology.')
        return None

    root_section = apical_root_sections[0]

    if method == 'apical_trunk':
        if apical_point is None:
            apical_secid = apical_point_section_segment(morph)[0]
        else:
            apical_secid = point_to_section_segment(morph, apical_point)[0] - 1

        apical_points = np.vstack(
            [section.points
             for section in morph.sections[apical_secid].iter(IterType.upstream)]
        )
        _points = []
        for section in root_section.iter():
            _points.append(section.points)
            if section.id == apical_secid:
                break
        apical_points = np.vstack(_points)
    elif method == 'first_segment':
        apical_points = root_section.points[:2]
    else:  # method == 'whole_apical'
        apical_points = np.vstack([section.points for section in root_section.iter()])

    if apical_points is None:
        logger.info('We did not find an apical point to align the morphology')
        return np.eye(3)
    if len(apical_points) < 2:
        logger.info('We did not find enough points in the apical to align the morphology')
        return np.eye(3)

    if len(apical_points) > 2:
        # compute covariance matrix
        cov = apical_points.T.dot(apical_points) / (len(apical_points) - 1)

        # compute principal direction
        principal_direction = svd(cov)[2][0]
        principal_direction *= np.sign(apical_points.dot(principal_direction).sum())
        principal_direction /= np.linalg.norm(principal_direction)
    else:
        principal_direction = apical_points[1] - apical_points[0]
        principal_direction /= np.linalg.norm(principal_direction)

    direction /= np.linalg.norm(direction)

    rotation_matrix = rotation_matrix_from_vectors(principal_direction, direction)

    # rotate the morphology in place
    rotate(morph, rotation_matrix)

    return rotation_matrix
