from pathlib import Path
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_equal

import morph_tool.nrnhines as test_module

DATA = Path(__file__).parent / 'data'
SIMPLE = DATA / 'simple2.asc'


def test_interpolate_compartments():
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [1, 2, 0],
                       [2, 2, 0]])
    paths = test_module._compartment_paths(points, 3)
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
    mapping = test_module.NeuroM_section_to_NRN_compartment_paths(SIMPLE)

    assert_array_almost_equal(mapping[1], [[[0., 0., 0.],
                                            [0., 5., 0.]]])
    assert_array_almost_equal(mapping[2], [[[0., 5., 0.],
                                            [-5., 5., 0.],
                                            [-6., 5., 0.]]])


def test_point_to_section_end():
    cell = test_module.get_NRN_cell(SIMPLE)
    assert_equal(test_module.point_to_section_end(cell.icell.all, [-8, 10, 0]),
                 6)

    assert_equal(test_module.point_to_section_end(cell.icell.all, [-8, 10, 2]),
                 None)

    assert_equal(test_module.point_to_section_end(cell.icell.all, [-8, 10, 2],
                                                  atol=10),
                 3)
