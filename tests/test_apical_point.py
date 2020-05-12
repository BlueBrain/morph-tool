import os


import neurom
from morphio import Morphology
from neurom.core.dataformat import COLS
from nose.tools import eq_
from morph_tool import apical_point_position, apical_point_section_segment


DATA = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

def test_get_apical_point():
    morph = Morphology(os.path.join(DATA, 'apical_test.swc'))
    point = apical_point_position(morph)
    eq_(point[COLS.Y], 25.0)

    # if we only take the very top, we only get one branch, whos common parents
    # is just itself
    point = apical_point_position(morph, tuft_percent=2)
    eq_(point[COLS.Y], 30.0)

    #try w/ a h5v2: this was converted using morphologyConverter
    morph = Morphology(os.path.join(DATA, 'apical_test.h5'))
    point = apical_point_position(morph)
    eq_(point[COLS.Y], 25.0)


def test__find_apical_section_segment():
    neuron = Morphology(os.path.join(DATA, 'apical_test.swc'))
    section, segment = apical_point_section_segment(neuron)
    eq_(section, 1)
    eq_(segment, 1)
