from pathlib import Path
import numpy as np
from numpy.testing import assert_array_almost_equal

import morph_tool.nrnhines as tested

DATA = Path(__file__).parent / 'data'
SIMPLE = DATA / 'simple2.asc'


def test_interpolate_compartments():
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [1, 2, 0],
                       [2, 2, 0]])
    paths = tested._compartment_paths(points, 3)
    assert len(paths) == 3

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

    assert_array_almost_equal(mapping[0], [[[0., 0., 0.],
                                            [0., 5., 0.]]])
    assert_array_almost_equal(mapping[1], [[[0., 5., 0.],
                                            [-5., 5., 0.],
                                            [-6., 5., 0.]]])


def test_point_to_section_end():
    cell = tested.get_NRN_cell(SIMPLE)
    assert (tested.point_to_section_end(cell.icell.all, [-8, 10, 0]) ==
                 6)

    assert (tested.point_to_section_end(cell.icell.all, [-8, 10, 2]) ==
                 None)

    assert (tested.point_to_section_end(cell.icell.all, [-8, 10, 2],
                                             atol=10) ==
                 3)


def test_NeuroM_section_to_NRN_section():
    # default use case
    mapping = tested.NeuroM_section_to_NRN_section(SIMPLE)
    assert mapping == {5: 0, 6: 1, 7: 2, 0: 3, 1: 4, 2: 5, 3: 6, 4: 7}

    # return reverse mapping
    mapping = tested.NeuroM_section_to_NRN_section(SIMPLE, reverse=True)
    assert mapping == {0: 5, 1: 6, 2: 7, 3: 0, 4: 1, 5: 2, 6: 3, 7: 4}

    # mapping with specific neurite type
    mapping = tested.NeuroM_section_to_NRN_section(SIMPLE, neurite_type="basal_dendrite")
    assert mapping == {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}

    # check we have a None for a zero len section
    mapping = tested.NeuroM_section_to_NRN_section(DATA / 'simple2_zero_len.asc')
    assert mapping == {5: 0, 6: 1, 7: 2, 0: 3, 1: None, 2: 4, 3: 5, 4: 6}

    # check if we have a soma, the sections id are correct
    mapping = tested.NeuroM_section_to_NRN_section(DATA / 'neuron.asc')
    assert mapping[79] == 1
    assert mapping[80] == 2


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
    assert result == [6] * 4
