"""Module to detect the terminal point of the main axon/"""
import logging
import numpy as np

from neurom import COLS
from neurom import NeuriteType, iter_sections
from morphio import SectionType, IterType

L = logging.getLogger(__name__)


def axonal_subtype(neuron, axonal_section=None, direction=None, bbox=None, ignore_axis=2):
    """Return a dict of extended axonal subtypes (main, collaterals)."""

    extended_types = dict()
    if axonal_section is None:
        axonal_section = neuron.sections[
            axon_point_section(neuron, direction=direction, ignore_axis=ignore_axis)
        ]
        for section in iter_sections(neuron):
            if section.type == SectionType.axon:
                extended_types[section.id] = "collateral"

        for section in axonal_section.iupstream():
            extended_types[section.id] = "main"
    else:

        for section in iter_sections(neuron):
            if section.type == SectionType.axon:
                extended_types[section.id] = "axon"
            if section.type == SectionType.basal_dendrite:
                extended_types[section.id] = "basal"
            if section.type == SectionType.apical_dendrite:
                extended_types[section.id] = "apical"

    return extended_types


def plot_axonal_subtypes(ax, nrn, subtype_map=None,
                         neurite_type=NeuriteType.all,
                         plane='xy',
                         soma_outline=True,
                         diameter_scale=None, linewidth=None,
                         color=None, alpha=None, realistic_diameters=False):
    """Same as neurom.view.plot_neuron, but plots subtypes of apicals.

    Args:
        subtype_map (dict): result of axonal_subtypes, if None, default axonal_subtype is called
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
        subtype_map = axonal_subtype(nrn)

    _plot_map = {
        "main": SectionType.custom5,
        "collateral": SectionType.custom6,
    }
    view.TREE_COLOR.update({
        NeuriteType.custom5: 'blue',
        NeuriteType.custom6: 'green',
    })

    for secid, subtype in subtype_map.items():
        nrn.sections[secid].morphio_section.type = _plot_map[subtype]

    view.plot_neuron(ax, nrn,
                     neurite_type=neurite_type,
                     plane=plane,
                     soma_outline=soma_outline,
                     diameter_scale=diameter_scale, linewidth=linewidth,
                     color=color, alpha=alpha, realistic_diameters=realistic_diameters)


def _get_angle(section, direction, axis):
    """Get angle between section endpoints and direction."""
    _diff = np.diff(section.points[:, axis], axis=0)
    angles = np.arccos(np.dot(direction[axis], _diff.T) / np.linalg.norm(_diff, axis=1))
    return list(angles[~np.isnan(angles)])  # in case we have duplicate points


def axon_point_section(morph, direction=None, bbox=None, ignore_axis=2):
    """Estimate axon point as the terminal point of the main axon.

    This point is defined as the point for which the sum of the angles with the direction vector
    is minimal. This corresponds to a main axon being a branch that follows best a straight line
    in the given direction. More constraints can be added with bbox argument, if incorect axon is
    detected.

    Args:
        morph (morphio.Morphology): morphology
        direction (ndarray): estimated direction of main axon, if None [0, -1, 0] is used
        bbox (dict): bbox of the form {'x': [-100, 100], 'y': [-500, -400]}
        ignore_axis (int): axis to ignore in the angle computation, if None, 3d directions are used

    Returns:
        MorphIO section ID for which the last point is the axon point
    """
    axis = [0, 1, 2]
    if ignore_axis is not None:
        axis.pop(ignore_axis)

    if direction is None:
        direction = np.array([0.0, -1.0, 0.0])
    else:
        direction = np.array(direction) / np.linalg.norm(direction)

    qualities = []
    ids = []
    for section in morph.iter():
        if section.type == SectionType.axon and not section.children:
            _angles = []
            for _section in section.iter(IterType.upstream):
                _angles += _get_angle(_section, direction, axis)
            qualities.append(np.mean(_angles))
            ids.append(section.id)

    if bbox is not None:
        qualities = np.array(qualities)
        positions = np.array([morph.sections[i].points[-1] for i in ids])
        ids = np.array(ids)
        _convert = {"x": COLS.X, "y": COLS.Y, "z": COLS.Z}
        for axis, bounds in bbox.items():
            bbox_filter = np.where(
                (positions[:, _convert[axis]] > bounds[0])
                & (positions[:, _convert[axis]] < bounds[1])
            )
            qualities = qualities[bbox_filter]
            positions = positions[bbox_filter]
            ids = ids[bbox_filter]

    if len(qualities) > 0:
        return ids[np.argmin(qualities)]
    else:
        L.warning("Could not find axon point in bounding box %s, we return None", bbox)
        return None
