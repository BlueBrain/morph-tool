import os
from nose.tools import eq_, assert_raises
from morphio import Morphology

from morph_tool import spatial


DATA = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

def test_point_to_section_segment():
    neuron = Morphology(os.path.join(DATA, 'apical_test.h5'))

    section, segment = spatial.point_to_section_segment(neuron, [0., 25., 0.])
    eq_(section, 2)
    eq_(segment, 1)

    assert_raises(ValueError, spatial.point_to_section_segment, neuron, [24, 0, 0])
