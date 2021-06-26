"""Module to retrieve the position of the apical point"""
import logging
import numpy as np

from neurom import COLS
from neurom import NeuriteType, iter_sections
from morphio import SectionType, IterType

from morph_tool.spatial import point_to_section_segment


L = logging.getLogger(__name__)

X, Y, Z = 0, 1, 2


def apical_subtype(neuron, apical_section=None, tuft_percent=20):
    """Return a dict of extended apical subtypes (trunk, oblique, tuft)."""
    extended_types = dict()

    apical_section = neuron.sections[
        apical_point_section_segment(neuron, tuft_percent=tuft_percent)[0]
    ]
    if apical_section is not None:
        for section in iter_sections(neuron):
            if section.type == SectionType.apical_dendrite:
                extended_types[section.id] = "oblique"

        for section in apical_section.ipreorder():
            extended_types[section.id] = "tuft"

        for section in apical_section.iupstream():
            extended_types[section.id] = "trunk"
    else:
        for section in iter_sections(neuron):
            if section.type == SectionType.apical_dendrite:
                extended_types[section.id] = "apical"

    for section in iter_sections(neuron):
        if section.type == SectionType.basal_dendrite:
            extended_types[section.id] = "basal"
        if section.type == SectionType.axon:
            extended_types[section.id] = "axon"

    return extended_types


def plot_apical_subtypes(ax, nrn, subtype_map=None,
                         neurite_type=NeuriteType.all,
                         plane='xy',
                         soma_outline=True,
                         diameter_scale=None, linewidth=None,
                         color=None, alpha=None, realistic_diameters=False):
    """Same as neurom.view.plot_neuron, but plots subtypes of apicals.

    Args:
        subtype_map (dict): result of apical_subtypes, if None, default apical_subtype is called
        see neurom.view.plot_neuron for all other parameters
    """

    from neurom.view import view

    if diameter_scale is None:
        diameter_scale = view._DIAMETER_SCALE
    if linewidth is None:
        linewidth = view._LINEWIDTH
    if alpha is None:
        alpha = view._ALPHA
    if subtype_map is None:
        subtype_map = apical_subtype(nrn)

    _plot_map = {
        "trunk": SectionType.custom5,
        "tuft": SectionType.custom6,
        "oblique": SectionType.custom7,
    }
    view.TREE_COLOR.update({
        NeuriteType.custom5: 'purple',
        NeuriteType.custom6: 'green',
        NeuriteType.custom7: 'orange',
    })

    for secid, subtype in subtype_map.items():
        nrn.sections[secid].morphio_section.type = _plot_map[subtype]

    view.plot_neuron(ax, nrn,
                     neurite_type=neurite_type,
                     plane=plane,
                     soma_outline=soma_outline,
                     diameter_scale=diameter_scale, linewidth=linewidth,
                     color=color, alpha=alpha, realistic_diameters=realistic_diameters)


def apical_point_section_segment(neuron, tuft_percent=20):
    """find the apical point's section and segment

    Args:
        neuron (morphio.Morphology): a morphology
        tuft_percent: percentage of the 'height' of the apical dendrite that
        would enclose the tuft, only leaves in this volume are considered as
        endpoints

    Returns:
        Tuple: (NeuroM/MorphIO section ID, point ID) of the apical point. Since NeuroM v2, section
        ids of NeuroM and MorphIO are the same excluding soma.
    """
    point = apical_point_position(neuron, tuft_percent=tuft_percent)

    if point is None:
        L.warning("Could not find apical point")
        return None, None

    section, segment = point_to_section_segment(neuron, point)
    return section, segment


def apical_point_position(neuron, tuft_percent=20):
    """Attempt to find the apical point in 'tufted' neurons

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
    """
    apical = [root for root in neuron.root_sections if SectionType.apical_dendrite == root.type]

    if not len(apical):
        return None
    assert len(apical) >= 1, "Too many apical dendrites"
    apical = apical[0]

    MIN, MAX = 0, 1

    points = np.vstack([section.points for section in apical.iter()])
    bounding_box = np.array([np.min(points[:, 0:3], axis=0), np.max(points[:, 0:3], axis=0)])
    if neuron.soma.center[Y] < apical.points[:, Y].mean():
        y_max = bounding_box[MAX][Y]
        y_min = (1 - tuft_percent / 100.0) * y_max
    else:
        y_min = bounding_box[MIN][Y]
        y_max = (1 - tuft_percent / 100.0) * y_min

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
