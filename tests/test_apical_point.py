import os


import neurom
from neurom.core.dataformat import COLS
from nose.tools import eq_
from morph_tool import get_apical_point


DATA_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

def test_get_apical_point():
    morph_path = os.path.join(DATA_PATH, 'apical_test.swc')
    morph = neurom.load_neuron(morph_path)
    point = get_apical_point(morph)
    eq_(point[COLS.Y], 25.0)

    # if we only take the very top, we only get one branch, whos common parents
    # is just itself
    point = get_apical_point(morph, tuft_percent=2)
    eq_(point[COLS.Y], 30.0)

    #try w/ a h5v2: this was converted using morphologyConverter
    morph_path = os.path.join(DATA_PATH, 'apical_test.h5')
    morph = neurom.load_neuron(morph_path)
    point = get_apical_point(morph)
    eq_(point[COLS.Y], 25.0)
