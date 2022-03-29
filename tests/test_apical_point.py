from pathlib import Path
from morphio import Morphology
from neurom.core.dataformat import COLS
from morph_tool import apical_point_position, apical_point_section_segment


DATA = Path(__file__).parent / 'data'

def test_get_apical_point():
    morph = Morphology(DATA / 'apical_test.swc')
    point = apical_point_position(morph)
    assert point[COLS.Y] == 25.0

    # if we only take the very top, we only get one branch, whose common parents
    # is just itself
    point = apical_point_position(morph, tuft_percent=2)
    assert point[COLS.Y] == 30.0

    #try w/ a h5v2: this was converted using morphologyConverter
    morph = Morphology(DATA / 'apical_test.h5')
    point = apical_point_position(morph)
    assert point[COLS.Y] == 25.0


def test__find_apical_section_segment():
    neuron = Morphology(DATA / 'apical_test.swc')
    section, segment = apical_point_section_segment(neuron)
    assert section == 1
    assert segment == 1
