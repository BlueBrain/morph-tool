from pathlib import Path
from morphio import Morphology, SectionType
from neurom import load_morphology, iter_sections
from morph_tool import sub_sectiontypes as tested


DATA = Path(__file__).absolute().parent / "data"


def test_apical_subtypes():
    morph = load_morphology(DATA / "neuron.asc")
    morph = tested.set_apical_subtypes(morph)
    types = [section.type for section in iter_sections(morph)]
    assert types[40] == SectionType.custom5
    assert types[47] == SectionType.custom6
    assert types[57] == SectionType.custom7

    morph = tested.unset_apical_subtypes(morph)
    types = [section.type for section in iter_sections(morph)]
    assert types[40] == SectionType.apical_dendrite
    assert types[47] == SectionType.apical_dendrite
    assert types[57] == SectionType.apical_dendrite


def test_axonal_subtypes():
    morph = load_morphology(DATA / "neuron.asc")
    morph = tested.set_axon_subtypes(morph)

    types = [section.type for section in iter_sections(morph)]
    assert types[79] == SectionType.custom5
    assert types[84] == SectionType.custom6

    morph = tested.unset_axon_subtypes(morph)
    types = [section.type for section in iter_sections(morph)]
    assert types[79] == SectionType.axon
    assert types[84] == SectionType.axon
