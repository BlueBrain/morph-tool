import os
import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_raises, assert_equal

from morph_tool.graft import graft_axon, find_axon
from morph_tool import MorphToolException

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
                       [[ 0. ,  0. ,  0. ],
                        [ 1. ,  2. ,  0. ],
                        [ 1. ,  2. ,  0. ],
                        [ 3. ,  4.5,  0. ],
                        [ 1. ,  2. ,  0. ],
                        [-5. , -4. ,  0. ],
                        [-5. , -4. ,  1. ]])


def test_self_graft():
    '''Grafting a neuron with its own neuron'''
    filename = os.path.join(_path, 'neuron.asc')
    new_axon = find_axon(ImmutMorphology(filename))

    neuron = Morphology(filename)
    graft_axon(neuron, new_axon)

    expected = Morphology(filename)
    assert_equal(expected, neuron)
