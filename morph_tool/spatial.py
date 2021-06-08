'''Module for spatial related functions'''
import numpy as np

from neurom import COLS

X, Y, Z = 0, 1, 2


def point_to_section_segment(neuron, point):
    '''Find section and segment that matches the point.

    Only the first point found with the *exact* same coordinates as the point argument is considered

    Args:
        neuron (morphio.Morphology): neuron object
        point (point): value of the point to find in the h5 file

    Returns:
        Tuple: (NeuroM/MorphIO section ID, point ID) of the point the matches the input coordinates.
        Since NeuroM v2, section ids of NeuroM and MorphIO are the same excluding soma.
    '''

    for section in neuron.iter():
        points = section.points
        offset = np.where((points[:, X] == point[COLS.X]) &
                          (points[:, Y] == point[COLS.Y]) &
                          (points[:, Z] == point[COLS.Z]))
        if offset[0].size:
            return section.id, offset[0][0]

    raise ValueError('Cannot find point in morphology that matches: {}'.format(point))
