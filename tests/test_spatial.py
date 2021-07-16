import os
from nose.tools import eq_, assert_raises
from pathlib import Path
from morphio import Morphology

from morph_tool import spatial


DATA = Path(__file__).absolute().parent / "data"


def test_point_to_section_segment():
    neuron = Morphology(os.path.join(DATA, 'apical_test.h5'))

    section, segment = spatial.point_to_section_segment(neuron, [0., 25., 0.])
    eq_(section, 1)
    eq_(segment, 1)

    assert_raises(ValueError, spatial.point_to_section_segment, neuron, [24, 0, 0])

    section, segment = spatial.point_to_section_segment(
        neuron, [0., 25.1, 0.], rtol=1e-1, atol=1e-1
    )
    eq_(section, 1)
    eq_(segment, 1)

    section, segment = spatial.point_to_section_segment(neuron, [0., 25.0001, 0.])
    eq_(section, 1)
    eq_(segment, 1)
