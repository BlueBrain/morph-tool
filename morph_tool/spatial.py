'''Module for spatial related functions'''
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.linalg import svd
from scipy.optimize import minimize_scalar

from morph_tool.transform import rotate
from neurom import COLS
from morphio import SectionType

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


def align_pyr_morph(morph, direction=[0.0, 1.0, 0.0]):
    """In place alignment of a pyramidal morphology such that apical tree is along 'direction'.

    The algorithm is based on SVD decomposititon of the correlation matrix obtained from all points
    in apical neurites, giving the principal axis of apical dendrite.
    If not apicals are present, no rotation is applied.

    Args:
        morph (morphio.Morphology): morphology to align
        direction (ndarray): 3 vector with final morphology direction
    """
    apical_points = []
    for section in morph.iter():
        if section.type == SectionType.apical_dendrite:
            apical_points += list(section.points)

    if len(apical_points) == 0:
        return None

    apical_points = np.array(apical_points)

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

    # rotate the morphology in place
    rotate(morph, rot(minimize_scalar(cost).x))
