from pathlib import Path
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_equal

import morph_tool.nrnhines as tested
from multiprocessing.pool import Pool

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


def _to_be_isolated(morphology_path, point):
    """Convert a point position to NEURON section index and return cell name and id.

    Args:
        morphology_path (Path): path to morphology
        point (list/ndarray): 3d point
    """
    cell = tested.get_NRN_cell(morphology_path)
    return tested.point_to_section_end(cell.icell.all, point)


def _isolated(morph_data):
    """Serves as a decorator to isolate convert_point_to_isc, but work with multiprocessing."""
    return tested.isolate(_to_be_isolated)(*morph_data)


def test_isolate():
    n_workers = 4
    with tested.NestedPool(processes=n_workers) as pool:
        result = list(pool.imap_unordered(_isolated,
                                          [(SIMPLE, [-8, 10, 0]),
                                           (SIMPLE, [-8, 10, 0]),
                                           (SIMPLE, [-8, 10, 0]),
                                           (SIMPLE, [-8, 10, 0]), ]))
    assert_equal(result, [6] * 4)
