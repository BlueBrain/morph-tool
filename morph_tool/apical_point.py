'''Module to retrieve the position of the apical point'''
import neurom
from neurom import geom, COLS

X, Y, Z = 0, 1, 2


def get_apical_point(morph, tuft_percent=20):
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
        morph is neurom.fst._core.Neuron
        tuft_percent: percentage of the 'height' of the apical dendrite that
        would enclose the tuft, only leaves in this volume are considered as
        endpoints

    Returns:
        neurom.core.point if point is found, or None if it isn't
    '''
    apical = [neurite for neurite in morph.neurites
              if neurom.NeuriteType.apical_dendrite == neurite.type]

    if not len(apical):
        return None
    assert len(apical) >= 1, 'Too many apical dendrites'
    apical = apical[0]

    MIN, MAX = 0, 1
    if morph.soma.center[Y] < apical.points[:, Y].mean():
        y_max = geom.bounding_box(apical)[MAX][Y]
        y_min = (1 - tuft_percent / 100.) * y_max
    else:
        y_min = geom.bounding_box(apical)[MIN][Y]
        y_max = (1 - tuft_percent / 100.) * y_min

    common_parents = set(apical.iter_sections())
    for leaf in apical.root_node.ileaf():
        end_point = leaf.points[-1, COLS.Y]
        if y_min <= end_point <= y_max:
            parents = leaf.iupstream()
            common_parents &= set(parents)

    for parent_section in apical.iter_sections():
        if parent_section in common_parents:
            common_parents.remove(parent_section)
            if not common_parents:
                return parent_section.points[-1]

    return None
