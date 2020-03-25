'''Module to retrieve the position of the apical point'''
import logging
import numpy as np

from neurom import COLS

from morphio import SectionType, IterType

from morph_tool.spatial import point_to_section_segment


logger = logging.getLogger('morph_tool')

X, Y, Z = 0, 1, 2


def apical_point_section_segment(neuron):
    '''find the apical point's section and segment

    Args:
        neuron (morphio.Morphology): a morphology

    Returns:
        a tuple (MorphIO section ID, point ID) of the apical point
    '''
    point = apical_point_position(neuron)

    if point is None:
        logger.warning('Could not find apical point')
        return None, None

    section, segment = point_to_section_segment(neuron, point)
    section -= 1  # MorphIO ID = NeuroM ID - 1
    return section, segment


def apical_point_position(neuron, tuft_percent=20):
    '''Attempt to find the apical point in 'tufted' neurons

    The algorithm is a simplification of:
https://bbpcode.epfl.ch/source/xref/analysis/Pneumatk/pneumatk/__tools__/ \
    Tree/methods/get_apical_point_index.py

    Consider a neuron:

        |   /    | Tuft = 20%
        |--/     |
        |   /
        |--/
        |
    ----.-----

    All endpoints in the top 'tuft_percent' are found, then their common
    branch segment, furthest from the soma, is identified.

    Using the release from 2012 as the base, apical points were compared from
    the annotated versions, and using this algorithm.  Of 239 morphologies that
    were annotated, 48 differed in the apical point choice by more than 1um
    in the y component.  Many of the differences were morphologies that
    shoulnd't have had apical points in the first place (ie: weren't pyramidal
    cells, ex: C050398B-I4.)

    Args:
        neuron (morphio.Morphology): a neuron
        tuft_percent: percentage of the 'height' of the apical dendrite that
        would enclose the tuft, only leaves in this volume are considered as
        endpoints

    Returns:
        neurom.core.point if point is found, or None if it isn't
    '''
    apical = [root for root in neuron.root_sections
              if SectionType.apical_dendrite == root.type]

    if not len(apical):
        return None
    assert len(apical) >= 1, 'Too many apical dendrites'
    apical = apical[0]

    MIN, MAX = 0, 1

    points = np.vstack([section.points for section in apical.iter()])
    bounding_box = np.array([np.min(points[:, 0:3], axis=0), np.max(points[:, 0:3], axis=0)])
    if neuron.soma.center[Y] < apical.points[:, Y].mean():
        y_max = bounding_box[MAX][Y]
        y_min = (1 - tuft_percent / 100.) * y_max
    else:
        y_min = bounding_box[MIN][Y]
        y_max = (1 - tuft_percent / 100.) * y_min

    common_parents = set(section.id for section in apical.iter())
    for leaf in apical.iter():
        if leaf.children:
            continue
        end_point = leaf.points[-1, COLS.Y]
        if y_min <= end_point <= y_max:
            parents = leaf.iter(IterType.upstream)
            common_parents &= set(section.id for section in parents)

    for parent_section in apical.iter():
        if parent_section.id in common_parents:
            common_parents.remove(parent_section.id)
            if not common_parents:
                return parent_section.points[-1]

    return None
