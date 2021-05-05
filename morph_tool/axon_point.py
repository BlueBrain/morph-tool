"""Module to detect the terminal point of the main axon/"""
import logging
import numpy as np

from neurom import COLS
from morphio import SectionType, IterType

L = logging.getLogger(__name__)


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
