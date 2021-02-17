'''Module for spatial related functions'''
import logging
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.linalg import svd
from scipy.optimize import minimize_scalar

from neurom import COLS
from morphio import SectionType, IterType

from morph_tool.transform import rotate

logger = logging.getLogger('morph_tool')
X, Y, Z = 0, 1, 2


def point_to_section_segment(neuron, point):
    '''find section and segment that matches the point

    Only the first point found with the *exact* same coordinates as the point argument is considered

    Args:
        neuron (morphio.Morphology): neuron object
        point (point): value of the point to find in the h5 file

    Returns a tuple (NeuroM section ID, point ID) of the point the matches the input coordinates
    '''

    for section in neuron.iter():
        points = section.points
        offset = np.where((points[:, X] == point[COLS.X]) &
                          (points[:, Y] == point[COLS.Y]) &
                          (points[:, Z] == point[COLS.Z]))
        if offset[0].size:
            # Section ids start at 0 in MorphIO
            return (section.id + 1, offset[0][0])

    raise ValueError('Cannot find point in morphology that matches: {}'.format(point))


def align_pyr_morph(morph, direction=None, method='whole_apical', apical_point=None):
    """In-place alignment of a pyramidal morphology such that apical tree is along 'direction'.

    The base algorithm is based on SVD decomposition of the correlation matrix obtained from
    points in apical neurites, giving the principal axis of apical dendrite.
    Currently, three algorithms are implemented, differing in the choice of points:

    1) with method='whole_apical': All the points in the apical dendrite are used.
    2) with method='first_section': Only the points in the first section are used.
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
    if method not in ['whole_apical', 'apical_trunk', 'first_section']:
        raise NotImplementedError(f"Method {method} is not implementd")

    if direction is None:
        direction = [0.0, 1.0, 0.0]
    apical_points = None
    for root_section in morph.root_sections:
        if root_section.type == SectionType.apical_dendrite:
            if apical_points is not None:
                logger.info('The morphology has two apical dendrites, we will use the second.')

            if method == 'apical_trunk':
                from morph_tool.apical_point import apical_point_section_segment
                if apical_point is None:
                    apical_secid = apical_point_section_segment(morph)[0]
                else:
                    from morph_tool.spatial import point_to_section_segment

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
            if method == 'first_section':
                apical_points = root_section.points
            else:
                apical_points = np.vstack([section.points for section in root_section.iter()])

    if apical_points is None:
        logger.info('We did not find an apical point to align the morphology')
        return np.eye(3)

    # compute covariance matrix
    cov = apical_points.T.dot(apical_points) / (len(apical_points) - 1)

    # compute principal direction
    principal_direction = svd(cov)[2][0]
    principal_direction *= np.sign(apical_points.dot(principal_direction).sum())
    principal_direction /= np.linalg.norm(principal_direction)
    direction /= np.linalg.norm(direction)

    def rot(angle):
        """Rotation matrix in the plane formed by principal_direction and direction."""
        return R.from_rotvec(angle * np.cross(principal_direction, direction)).as_matrix()

    def cost(angle):
        """Cost to find the rotation angle."""
        return np.linalg.norm(rot(angle).dot(principal_direction) - direction)

    rotation_matrix = rot(minimize_scalar(cost).x)

    # rotate the morphology in place
    rotate(morph, rotation_matrix)

    return rotation_matrix
