"""Module for spatial related functions."""
import numpy as np

from neurom import COLS


def point_to_section_segment(neuron, point, rtol=1e-05, atol=1e-08):
    """Find section and segment that matches the point.

    Only the first point found with the *exact* same coordinates as the point argument is considered

    Args:
        neuron (morphio.Morphology): neuron object
        point (point): value of the point to find in the h5 file
        rtol, atol (floats): precision of np.isclose

    Returns:
        Tuple: (NeuroM/MorphIO section ID, point ID) of the point the matches the input coordinates.
        Since NeuroM v2, section ids of NeuroM and MorphIO are the same excluding soma.
    """
    for section in neuron.iter():
        points = section.points
        offset = np.where(
            np.isclose(points[:, COLS.XYZ], point[COLS.XYZ], rtol=rtol, atol=atol).all(axis=1)
        )
        if offset[0].size:
            return section.id, offset[0][0]

    raise ValueError(f'Cannot find point in morphology that matches: {point}')
