from pathlib import Path
from morphio import Morphology
from morph_tool import axon_point_section


DATA = Path(__file__).absolute().parent / "data"


def test_axon_point_section():
    morph = Morphology(DATA / "neuron.asc")

    assert axon_point_section(morph) == 98
    assert axon_point_section(morph, direction=[1.0, 0, 0]) == 140
    assert axon_point_section(morph, bbox={"x": [-10, 10]}) == 136
    assert axon_point_section(morph, bbox={"x": [-100, 100], "y": [-100, 100]}) == 142
