from pathlib import Path
from pkg_resources import get_distribution, parse_version

from morph_tool import diff
from morphio import PointLevel, SectionType, Warning, set_ignored_warning
from morphio.mut import Morphology

DATA = Path(__file__).parent / 'data'


def test_equality():
    set_ignored_warning([Warning.wrong_duplicate, Warning.only_child], True)
    filename = DATA / 'simple2.asc'
    neuron_ref = Morphology(filename)
    assert not diff(neuron_ref, neuron_ref)
    assert not diff(neuron_ref, filename)
    assert not diff(neuron_ref, neuron_ref.as_immutable())

    def mundane_section(neuron):
        '''Not a root section, not a leaf section'''
        return neuron.root_sections[0].children[0]

    # Report different section types
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).type = SectionType.apical_dendrite
    result = diff(neuron_ref, a)
    assert result
    assert result.info == 'Section(id=1, points=[(0 5 0),..., (-6 5 0)]) and Section(id=1, points=[(0 5 0),..., (-6 5 0)]) have different section types'

    # Report different section points
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).points = [[0,0,0], [0,0,0], [0,0,1]]
    result = diff(neuron_ref, a)
    assert result
    assert (result.info ==
                 '\n'.join(['Attributes Section.points of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 0 0),..., (0 0 1)])',
                            'have the same shape but different values',
                            'Vector points differs at index 0: [0. 5. 0.] != [0. 0. 0.]']))

    # Report different section diameters
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).diameters = [0,0,0]
    result = diff(neuron_ref, a)
    assert result
    assert (result.info ==
                 '\n'.join(['Attributes Section.diameters of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'have the same shape but different values',
                            'Vector diameters differs at index 0: 3.0 != 0.0']))

    # Report different section perimeters
    a = Morphology(DATA / 'simple2.asc')
    for section in a.iter():
        section.perimeters = [1] * len(section.points)
    result = diff(neuron_ref, a)
    assert result
    assert (result.info ==
                 '\n'.join(['Attributes Section.perimeters of:',
                            'Section(id=0, points=[(0 0 0),..., (0 5 0)])',
                            'Section(id=0, points=[(0 0 0),..., (0 5 0)])',
                            'have different shapes: (0,) vs (2,)']))


    # Report different number of children
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).append_section(PointLevel([[-6, 5, 0], [4, 5, 6]], [2, 3]))
    result = diff(neuron_ref, a)
    assert result
    assert result.info == 'Section(id=1, points=[(0 5 0),..., (-6 5 0)]) and Section(id=1, points=[(0 5 0),..., (-6 5 0)]) have a different number of children'

    # Report different number of root sections
    a = Morphology(DATA / 'simple2.asc')
    a.delete_section(a.root_sections[0])
    result = diff(neuron_ref, a)
    assert result
    assert result.info == 'Both morphologies have a different number of root sections'

    # Report only different number of root sections even though diameters are also different
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).diameters = [0,0,0]
    a.delete_section(a.root_sections[1])
    result = diff(neuron_ref, a, all_diffs=True)
    assert result
    assert result.info == 'Both morphologies have a different number of root sections'

    # Report only different section points but not the different section diameters
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).points = [[0,0,0], [0,0,0], [0,0,1]]
    mundane_section(a).diameters = [0,0,0]
    result = diff(neuron_ref, a, all_diffs=False)
    assert result
    assert (result.info ==
                 '\n'.join(['Attributes Section.points of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 0 0),..., (0 0 1)])',
                            'have the same shape but different values',
                            'Vector points differs at index 0: [0. 5. 0.] != [0. 0. 0.]']))

    # Report both different section points and different section diameters
    a = Morphology(DATA / 'simple2.asc')
    mundane_section(a).points = [[0,0,0], [0,0,0], [0,0,1]]
    mundane_section(a).diameters = [0,0,0]
    result = diff(neuron_ref, a, all_diffs=True)
    assert result
    assert (result.info ==
                 '\n'.join(['Attributes Section.points of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 0 0),..., (0 0 1)])',
                            'have the same shape but different values',
                            'Vector points differs at index 0: [0. 5. 0.] != [0. 0. 0.]',
                            '',
                            'Attributes Section.diameters of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 0 0),..., (0 0 1)])',
                            'have the same shape but different values',
                            'Vector diameters differs at index 0: 3.0 != 0.0']))

    # Check result for MorphIO < 3
    result = diff(DATA / 'single_child.asc', DATA / 'not_single_child.asc')
    if parse_version(get_distribution("morphio").version) < parse_version("3"):
        assert not result
    else:
        assert result
    set_ignored_warning([Warning.wrong_duplicate, Warning.only_child], False)
