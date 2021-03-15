import os

from morphio import Morphology
from nose.tools import eq_
from morph_tool import axon_point_section


DATA = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")


def test_get_apical_point():
    morph = Morphology(os.path.join(DATA, "neuron.asc"))

    eq_(axon_point_section(morph), 100)
    eq_(axon_point_section(morph, direction=[1.0, 0, 0]), 142)
    eq_(axon_point_section(morph, bbox={"x": [-10, 10]}), 136)
    eq_(axon_point_section(morph, bbox={"x": [-100, 100], "y": [-100, 100]}), 100)
