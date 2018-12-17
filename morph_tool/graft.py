'''What we call grafting is the process of taking one piece of a neuron and
putting it on another neuron.'''
from morphio import SectionType

from morph_tool import MorphToolException
from morph_tool.transform import translate


def find_axon(neuron):
    '''Find the root section which is an axon

    raises if number of axon found is not 1'''
    axons = [root for root in neuron.root_sections
             if root.type == SectionType.axon]

    if not axons:
        raise MorphToolException('No axon found!')

    if len(axons) > 1:
        raise MorphToolException('Too many axons found!')

    return axons[0]


def graft_axon(neuron, axon):
    '''Graft an axon

    This is a very simple implementation.
    No check are performed and the new axon is simply
    translated at the place where the old axon was.
    '''

    old_axon = find_axon(neuron)
    axon_start = old_axon.points[0]

    neuron.delete_section(old_axon)
    new_axon = neuron.append_root_section(axon, recursive=True)

    translate(new_axon, axon_start - new_axon.points[0])
