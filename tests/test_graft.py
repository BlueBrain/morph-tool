import os
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_raises, assert_equal, assert_almost_equal

from morph_tool.graft import graft_axon, find_axon, _dendrites_mean_direction, _axon_dendrites_angle, _random_direction, _rotation_matrix
from morph_tool import MorphToolException, NoAxonException

from morphio.mut import Morphology
from morphio import Morphology as ImmutMorphology, Option

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
SIMPLE = Morphology(os.path.join(_path, 'simple.swc'))

def test_find_axon():
    axon = find_axon(SIMPLE)
    assert_array_equal(axon.points,
                       [[0, 0, 0], [0, -4, 0]])


    with assert_raises(MorphToolException) as obj:
        find_axon(Morphology(os.path.join(_path, 'no_axon.swc')))
    assert_equal(str(obj.exception), 'No axon found!')

    with assert_raises(MorphToolException) as obj:
        find_axon(Morphology(os.path.join(_path, 'two_axons.asc')))
    assert_equal(str(obj.exception), 'Too many axons found!')


def test_graft():
    m = Morphology(os.path.join(_path, 'simple2.swc'))
    new_axon = find_axon(m)

    neuron = Morphology(os.path.join(_path, 'simple.swc'))
    graft_axon(neuron, new_axon)

    grafted_axon = find_axon(neuron)
    points = np.vstack([section.points for section in grafted_axon.iter()])

    assert_array_equal(points,
                       np.array([[ 0.0000000e+00,  0.0000000e+00,  0.0000000e+00],
                                 [-2.4447479e-07, -2.2360680e+00,  0.0000000e+00],
                                 [-2.4447479e-07, -2.2360680e+00,  0.0000000e+00],
                                 [-6.7082095e-01, -5.3665628e+00,  0.0000000e+00],
                                 [-2.4447479e-07, -2.2360680e+00,  0.0000000e+00],
                                 [ 2.6832821e+00,  5.8137765e+00,  0.0000000e+00],
                                 [ 2.6832821e+00,  5.8137765e+00,  1.0000000e+00]],
                                dtype=np.float32))



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
                              [[5.110419 , 5.486378 , 4.9647303],
                               [5.9937673, 9.377404 , 4.6825743]])

def test_self_graft():
    '''Grafting a neuron with its own neuron'''
    filename = os.path.join(_path, 'neuron.asc')
    new_axon = find_axon(ImmutMorphology(filename))

    neuron = Morphology(filename)
    graft_axon(neuron, new_axon)

    expected = Morphology(filename)
    assert_equal(expected, neuron)


def test_rotation_matrix():
    r = _rotation_matrix([1, 0, 0], [0, 1, 0])
    assert_array_almost_equal(r.dot([1,0,0]),
                              [0, 1, 0])

    u = [1, 2, 4]
    v = [3, 0, 8]
    u /= np.linalg.norm(u)
    v /= np.linalg.norm(v)

    r = _rotation_matrix(u, v)
    assert_array_almost_equal(r.dot(u), v)
