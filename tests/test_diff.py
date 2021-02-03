from os.path import dirname
from os.path import join as joinp

from click.testing import CliRunner
from morph_tool import diff
from morphio import PointLevel, SectionType, Warning, set_ignored_warning
from morphio.mut import Morphology
from nose.tools import assert_equal, ok_

DATA = joinp(dirname(__file__), 'data')


def test_equality():
    set_ignored_warning([Warning.wrong_duplicate, Warning.only_child], True)
    filename = joinp(DATA, 'simple2.asc')
    neuron_ref = Morphology(filename)
    ok_(not diff(neuron_ref, neuron_ref))
    ok_(not diff(neuron_ref, filename))
    ok_(not diff(neuron_ref, neuron_ref.as_immutable()))

    def mundane_section(neuron):
        '''Not a root section, not a leaf section'''
        return neuron.root_sections[0].children[0]

    a = Morphology(joinp(DATA, 'simple2.asc'))
    mundane_section(a).type = SectionType.apical_dendrite
    result = diff(neuron_ref, a)
    ok_(result)
    assert_equal(result.info, 'Section(id=1, points=[(0 5 0),..., (-6 5 0)]) and Section(id=1, points=[(0 5 0),..., (-6 5 0)]) have different section types')

    a = Morphology(joinp(DATA, 'simple2.asc'))
    mundane_section(a).points = [[0,0,0], [0,0,0], [0,0,1]]
    result = diff(neuron_ref, a)
    ok_(result)
    assert_equal(result.info,
                 '\n'.join(['Attributes Section.points of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 0 0),..., (0 0 1)])',
                            'have the same shape but different values',
                            'Vector points differs at index 0: [0. 5. 0.] != [0. 0. 0.]']))

    a = Morphology(joinp(DATA, 'simple2.asc'))
    mundane_section(a).diameters = [0,0,0]
    result = diff(neuron_ref, a)
    ok_(result)
    assert_equal(result.info,
                 '\n'.join(['Attributes Section.diameters of:',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'Section(id=1, points=[(0 5 0),..., (-6 5 0)])',
                            'have the same shape but different values',
                            'Vector diameters differs at index 0: 3.0 != 0.0']))

    a = Morphology(joinp(DATA, 'simple2.asc'))
    for section in a.iter():
        section.perimeters = [1] * len(section.points)
    result = diff(neuron_ref, a)
    ok_(result)
    assert_equal(result.info,
                 '\n'.join(['Attributes Section.perimeters of:',
                            'Section(id=0, points=[(0 0 0),..., (0 5 0)])',
                            'Section(id=0, points=[(0 0 0),..., (0 5 0)])',
                            'have different shapes: (0,) vs (2,)']))


    a = Morphology(joinp(DATA, 'simple2.asc'))
    mundane_section(a).append_section(PointLevel([[-6, 5, 0], [4, 5, 6]], [2, 3]))
    result = diff(neuron_ref, a)
    ok_(result)
    assert_equal(result.info, 'Section(id=1, points=[(0 5 0),..., (-6 5 0)]) and Section(id=1, points=[(0 5 0),..., (-6 5 0)]) have different a different number of children')

    a = Morphology(joinp(DATA, 'simple2.asc'))
    a.delete_section(a.root_sections[0])
    result = diff(neuron_ref, a)
    ok_(result)
    assert_equal(result.info, 'Both morphologies have different a different number of root sections')

    result = diff(joinp(DATA, 'single_child.asc'), joinp(DATA, 'not_single_child.asc'))
    ok_(not result)
    set_ignored_warning([Warning.wrong_duplicate, Warning.only_child], False)
