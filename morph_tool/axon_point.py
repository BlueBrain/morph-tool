"""Module to detect the terminal point of the main axon/"""
import numpy as np

from morphio import SectionType, IterType
import logging

logger = logging.getLogger("morph_tool")


def get_axon_point(morph, direction=None, bbox=None):
    """Estimate axon point as the terminal point of the main axon.

    This point is defined as the point for which the sum of the angles with the direction vector
    is minimal. This corresponds to a main axon being a branch that follows best a straight line
    in the given direction. More constraints can be added with bbox argument, if incorect axon is
    detected.

    Args:
        morph (morphio.Morphology): morphology
        direction (ndarray): estimated direction of main axon, if None [0, -1, 0] is used
        bbox (dict): bbox of the form {'x': [-100, 100], 'y': [-500, -400]}

    Returns:
        MorphIO section ID for which the last point is the axon point
    """
    if direction is None:
        direction = [0.0, -1.0, 0.0]
    else:
        direction /= np.linalg.norm(direction)

    def _get_angle(section, direction):
        """Get angle between section endpoints and direction."""
        return np.arccos(
            np.dot(direction, section.points[-1] - section.points[0])
            / (np.linalg.norm(section.points[-1] - section.points[0]))
        )

    qualities, ids, positions = [], [], []
    for section in morph.iter():
        if section.type == SectionType.axon and not section.children:
            qualities.append(
                sum(_get_angle(section, direction) for section in section.iter(IterType.upstream))
            )
            ids.append(section.id)
            positions.append(section.points[-1])

    if bbox is not None:
        qualities, positions, ids = np.array(qualities), np.array(positions), np.array(ids)
        _convert = {"x": 0, "y": 1, "z": 2}
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
        logger.warning("Could not find axon point in bounding box %s, we return None", bbox)
        return None
