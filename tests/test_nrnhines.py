from pathlib import Path
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_equal
from nose.tools import assert_dict_equal

import morph_tool.nrnhines as tested

DATA = Path(__file__).parent / 'data'
SIMPLE = DATA / 'simple2.asc'


def test_interpolate_compartments():
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [1, 2, 0],
                       [2, 2, 0]])
    paths = tested._compartment_paths(points, 3)
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


def test_NeuroM_section_to_NRN_section():
    mapping = tested.NeuroM_section_to_NRN_section(SIMPLE)
    assert_dict_equal(mapping,
                      {6: 0, 7: 1, 8: 2, 1: 3, 2: 4, 3: 5, 4: 6, 5: 7})

    mapping = tested.NeuroM_section_to_NRN_section(DATA / 'real.asc')
    assert_dict_equal(mapping,
                      dict(zip(range(1, 661), range(1, 661))))

def test_NeuroM_section_to_NRN_compartment_paths():
    mapping = tested.NeuroM_section_to_NRN_compartment_paths(SIMPLE)

    assert_array_almost_equal(mapping[1], [[[0., 0., 0.],
                                            [0., 5., 0.]]])
    assert_array_almost_equal(mapping[2], [[[0., 5., 0.],
                                            [-5., 5., 0.],
                                            [-6., 5., 0.]]])


def test_point_to_section_end():
    cell = tested.get_NRN_cell(SIMPLE)
    assert_equal(tested.point_to_section_end(cell.icell.all, [-8, 10, 0]),
                 6)

    assert_equal(tested.point_to_section_end(cell.icell.all, [-8, 10, 2]),
                 None)

    assert_equal(tested.point_to_section_end(cell.icell.all, [-8, 10, 2],
                                                  atol=10),
                 3)
