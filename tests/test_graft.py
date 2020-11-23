import os
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_raises, assert_equal, assert_almost_equal, ok_

from morph_tool.graft import graft_axon, find_axon, _dendrites_mean_direction, _axon_dendrites_angle, _random_direction, _soma_mean_radius
import morph_tool.graft as graft
from morph_tool import MorphToolException, NoAxonException, diff

from morphio.mut import Morphology
from morphio import PointLevel, SomaType, SectionType, Morphology as ImmutMorphology, Option

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
SIMPLE = Morphology(os.path.join(_path, 'simple.swc'))


def test_find_axon():
    axon = find_axon(SIMPLE)
    assert_array_equal(axon.points,
                       [[0, 0, 0], [0, -4, 0]])

    # section is not an axon
    with assert_raises(MorphToolException):
        dendrite = SIMPLE.root_sections[0]
        find_axon(dendrite)

    # no axon found
    with assert_raises(MorphToolException) as obj:
        find_axon(Morphology(os.path.join(_path, 'no_axon.swc')))
    assert_equal(str(obj.exception), 'No axon found!')

    # first axon is chosen
    axon = find_axon(Morphology(os.path.join(_path, 'two_axons.asc')))
    assert_array_equal(axon.points, [[0, 0, 0], [0, -4, 0]])


def test_graft():
    m = Morphology(os.path.join(_path, 'simple2.swc'))
    new_axon = find_axon(m)

    neuron = Morphology(os.path.join(_path, 'simple.swc'))
    graft_axon(neuron, new_axon)

    grafted_axon = find_axon(neuron)
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
    direction = _random_direction(axis, 2.3)
    assert_almost_equal(np.arccos(direction.dot(axis) / np.linalg.norm(direction) / np.linalg.norm(axis)), 2.3)

    axis2 = np.array([0, 3, 1])
    direction = _random_direction(axis2, 2.3)
    assert_almost_equal(np.arccos(direction.dot(axis2) / np.linalg.norm(direction) / np.linalg.norm(axis2)), 2.3)


def test_dendrites_mean_direction():
    donor_neuron = Morphology(os.path.join(_path, 'simple3.asc'))
    assert_array_equal(_dendrites_mean_direction(donor_neuron),
                       np.array([0.6666667, 4.6666665, 3.       ],
                                dtype=np.float32))

    no_dendrite = Morphology()
    no_dendrite.append_root_section(PointLevel([[0,0,0], [1,1,1]], [1, 1]), SectionType.axon)
    with assert_raises(MorphToolException) as obj:
        _dendrites_mean_direction(no_dendrite)

def test_soma_mean_radius():
    m = Morphology()
    m.soma.points = [[0,0,0], [1,0,0], [1,1,0], [0,1,0]]
    assert_equal(_soma_mean_radius(m, [0.5, 0.5, 0]),
                 0.7071067811865476)


def test_axon_dendrites_angle():
    assert_almost_equal(_axon_dendrites_angle(SIMPLE), np.pi, places=5)


def test_graft_axon_on_synthesized_cell():
    np.random.seed(0)
    # donor neuron is empty
    assert_raises(NoAxonException, graft_axon, SIMPLE, Morphology())

    donor_neuron = Morphology(os.path.join(_path, 'simple3.asc'))
    synthesized_cell = Morphology(os.path.join(_path, 'synthesized_cell.asc'))
    graft_axon(synthesized_cell, donor_neuron)
    axon = find_axon(synthesized_cell)
    assert_array_almost_equal(axon.points,
                              [[5.110419, 5.486378, 4.9647303],
                               [5.110419, 1.486378, 4.9647303]])

def test_self_graft():
    '''Grafting a neuron with its own neuron'''
    filename = os.path.join(_path, 'neuron.asc')
    new_axon = find_axon(ImmutMorphology(filename))

    neuron = Morphology(filename)
    graft_axon(neuron, new_axon)

    expected = Morphology(filename)
    ok_(not diff(expected, neuron))


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
