import numpy as np
import os
from numpy.testing import assert_array_almost_equal, assert_equal

from morph_tool.nrnhines import NeuroM_section_to_NRN_compartment_paths, _compartment_paths

PATH = os.path.join(os.path.dirname(__file__), 'data')


def test_interpolate_compartments():
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [1, 2, 0],
                       [2, 2, 0]])
    paths = _compartment_paths(points, 3)
    assert_equal(len(paths), 3)

    assert_array_almost_equal(paths[0],
                              [[0., 0., 0.],
                               [1., 0., 0.],
                               [1., 0.33333333, 0.]])

    assert_array_almost_equal(paths[1],
                              [[1., 0.33333333, 0.],
                               [1., 1.66666667, 0.]])

    assert_array_almost_equal(paths[2],
                              [[1., 1.66666667, 0.],
                               [1., 2., 0.],
                               [2., 2., 0.]])


def test_NeuroM_section_to_NRN_compartment_paths():
    filename = os.path.join(PATH, 'simple2.asc')
    mapping = NeuroM_section_to_NRN_compartment_paths(filename)

    assert_array_almost_equal(mapping[1], [[[0., 0., 0.],
                                            [0., 5., 0.]]])
    assert_array_almost_equal(mapping[2], [[[0., 5., 0.],
                                            [-5., 5., 0.],
                                            [-6., 5., 0.]]])
