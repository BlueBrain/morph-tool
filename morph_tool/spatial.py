'''Module for spatial related functions'''
import numpy as np
from neurom import COLS

X, Y, Z = 0, 1, 2


def point_to_section_segment(neuron, point):
    '''find section and segment that matches the point

    Only the first point found with the *exact* same coordinates as the point argument is considered

    Args:
        neuron (morphio.Morphology): path to file to examine
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
