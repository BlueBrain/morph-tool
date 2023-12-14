from pathlib import Path
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_almost_equal
import pytest

from morph_tool import diff, graft
from morph_tool.exceptions import MorphToolException, NoAxonException

from morphio.mut import Morphology
from morphio import PointLevel, SectionType, Morphology as ImmutMorphology

DATA = Path(__file__).parent / 'data'
SIMPLE = Morphology(DATA / 'simple.swc')


def test_find_axon():
    axon = graft.find_axon(SIMPLE)
    assert_array_equal(axon.points,
                       [[0, 0, 0], [0, -4, 0]])

    # section is not an axon
    with pytest.raises(MorphToolException):
        dendrite = SIMPLE.root_sections[0]
        graft.find_axon(dendrite)

    # no axon found
    with pytest.raises(MorphToolException, match='No axon found!'):
        graft.find_axon(Morphology(DATA / 'no_axon.swc'))

    # first axon is chosen
    axon = graft.find_axon(Morphology(DATA / 'two_axons.asc'))
    assert_array_equal(axon.points, [[0, 0, 0], [0, -4, 0]])


def test_graft():
    for rng in [np.random, np.random.default_rng(0)]:
        m = Morphology(DATA / 'simple2.swc')
        new_axon = graft.find_axon(m)

        neuron = Morphology(DATA / 'simple.swc')
        graft.graft_axon(neuron, new_axon)

        grafted_axon = graft.find_axon(neuron)
        points = np.vstack([section.points for section in grafted_axon.iter()])

        assert_array_equal(points,
                           np.array([[ 0. ,  0. ,  0. ],
                                     [ 1. ,  2. ,  0. ],
                                     [ 1. ,  2. ,  0. ],
                                     [ 3. ,  4.5,  0. ],
                                     [ 1. ,  2. ,  0. ],
                                     [-5. , -4. ,  0. ],
                                     [-5. , -4. ,  1. ]], dtype=np.float32))


def test_random_direction():
    axis = np.array([0, 0, 2])
    direction = graft._random_direction(axis, 2.3)
    assert_almost_equal(np.arccos(direction.dot(axis) / np.linalg.norm(direction) / np.linalg.norm(axis)), 2.3)

    axis2 = np.array([0, 3, 1])
    direction = graft._random_direction(axis2, 2.3)
    assert_almost_equal(np.arccos(direction.dot(axis2) / np.linalg.norm(direction) / np.linalg.norm(axis2)), 2.3)


def test_dendrites_mean_direction():
    donor_neuron = Morphology(DATA / 'simple3.asc')
    assert_array_equal(graft._dendrites_mean_direction(donor_neuron),
                       np.array([0.6666667, 4.6666665, 3.       ],
                                dtype=np.float32))

    no_dendrite = Morphology()
    no_dendrite.append_root_section(PointLevel([[0,0,0], [1,1,1]], [1, 1]), SectionType.axon)
    with pytest.raises(MorphToolException):
        graft._dendrites_mean_direction(no_dendrite)

def test_soma_mean_radius():
    m = Morphology()
    m.soma.points = [[0,0,0], [1,0,0], [1,1,0], [0,1,0]]
    assert (graft._soma_mean_radius(m, [0.5, 0.5, 0]) == 0.7071067811865476)


def test_axon_dendrites_angle():
    assert_almost_equal(graft._axon_dendrites_angle(SIMPLE), np.pi, decimal=5)


def test_graft_axon_on_synthesized_cell():
    np.random.seed(0)
    # donor neuron is empty
    with pytest.raises(NoAxonException):
        graft.graft_axon(SIMPLE, Morphology())

    donor_neuron = Morphology(DATA / 'simple3.asc')
    synthesized_cell = Morphology(DATA / 'synthesized_cell.asc')
    graft.graft_axon(synthesized_cell, donor_neuron)
    axon = graft.find_axon(synthesized_cell)
    assert_array_almost_equal(axon.points,
                              [[5.1272364, 5.4825425, 4.9689593],
                               [5.1272364, 1.4825425, 4.9689593]])


def test_self_graft():
    '''Grafting a neuron with its own neuron'''
    filename = DATA / 'neuron.asc'
    new_axon = graft.find_axon(ImmutMorphology(filename))

    neuron = Morphology(filename)
    graft.graft_axon(neuron, new_axon)

    expected = Morphology(filename)
    assert not diff(expected, neuron)


def test_hotfix_h5_duplicate():
    class Section:
        def __init__(self, points):
            self.points = np.array(points)

    # normal case
    assert_array_equal(graft._section_initial_direction(Section([[1, 2, 3], [4, 5, 6]])),
                       [3, 3, 3])
    # duplicate case
    assert_array_equal(graft._section_initial_direction(Section([[1, 2, 3],
                                                                 [1, 2, 3],
                                                                 [4, 5, 6]])),
                       [3, 3, 3])
