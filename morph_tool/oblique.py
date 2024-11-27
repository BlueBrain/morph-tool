"""Tools to manipulate oblique branches."""
import numpy as np
from morphio import SectionType, PointLevel
from neurom.morphmath import interval_lengths
from neurom.core import Morphology
from neurom import iter_sections

from morph_tool.sub_sectiontypes import (
    APICAL_SUBTYPE_MAP,
    set_apical_subtypes,
    unset_apical_subtypes,
)


def get_trunk_obliques_ids(morphology):
    """Return the ids of trunk of obliques (first oblique section off the trunk).

    This function assume that set_subtypes was run to know the obliques.
    """
    trunk_oblique_ids = []
    for section in iter_sections(morphology):
        if (
            section.type == APICAL_SUBTYPE_MAP["oblique"]
            and section.parent.type == APICAL_SUBTYPE_MAP["trunk"]
        ):
            trunk_oblique_ids.append(section.id)
    return trunk_oblique_ids


def remove_obliques(morphology, trunk_oblique_ids=None):
    """Remove obliqeus from morphology, from a list of trunk ids, or all is None."""
    morphology = Morphology(morphology)
    if trunk_oblique_ids is None:
        trunk_oblique_ids = get_trunk_obliques_ids(morphology)

    for trunk_id in trunk_oblique_ids:
        section = morphology.section(trunk_id)
        morphology.delete_section(section, recursive=True)
    return morphology


def remove_tuft(morphology):
    """Remove obliqeus from morphology, from a list of trunk ids, or all is None."""
    morphology = Morphology(morphology)
    morphology = set_apical_subtypes(morphology)
    for sec in iter_sections(morphology):
        if sec.type == APICAL_SUBTYPE_MAP["tuft"]:
            section = morphology.section(sec.id)
            morphology.delete_section(section, recursive=False)
    return unset_apical_subtypes(morphology)


def replace_apical(morphology, diameter=4.0, length=800):
    """Replace apical tree with stub apical consisting of 2 points."""
    for root_section in morphology.root_sections:
        if root_section.type == SectionType.apical_dendrite:
            orig_point = root_section.points[0]
            morphology.delete_section(root_section, recursive=True)

    points = np.zeros((2, 3))
    points[:, 1] = np.linspace(0, length, 2)
    points += orig_point
    morphology.append_root_section(PointLevel(points, 2 * [diameter]), APICAL_SUBTYPE_MAP["trunk"])


def add_oblique(morphology, distance=100, length=200.0, diameter=2.0, direction=None):
    """Add a stub oblique to a trunk.

    Args:
        morphology (morphio.mut.Morphology): morph to add oblique with custom5 trunk
        distance (float): distance of oblique along trunk from trunk first point
        length (float): lenght of oblique
        diameter (float): diameter of oblique
        direction (-1/1/None): left (-1), right (1) direction along x, None will be a random choice
    """

    def _copy(section, section_base):
        """to recursively copy downstream from section_base to section"""
        for base_child in section_base.children:
            section.append_section(base_child)
        for child, base_child in zip(section.children, section_base.children):
            _copy(child, base_child)

    oblique_points = np.zeros((2, 3))
    if direction is None:
        direction = (-1.0) ** np.random.randint(2)
    oblique_points[:, 0] = direction * np.linspace(0, length, 2)
    oblique_diameters = 2 * [diameter]
    current_dist = 0
    for section in morphology.iter():
        if section.type == APICAL_SUBTYPE_MAP["trunk"]:
            dists = current_dist + np.cumsum(interval_lengths(section.points, prepend_zero=True))
            current_dist = dists[-1]
            if dists[0] < distance <= dists[-1]:
                if distance == dists[-1]:
                    distance -= 0.1
                    print("we avoid trifurcation by moving oblique back by 0.1 micron")

                bif_pt = 0.5 * (
                    section.points[dists < distance][-1] + section.points[dists > distance][0]
                )
                bif_diam = 0.5 * (
                    section.diameters[dists < distance][-1] + section.diameters[dists > distance][0]
                )

                prev_points = np.append(section.points[dists < distance], [bif_pt], axis=0)
                prev_diams = np.append(section.diameters[dists < distance], bif_diam)
                prev_pts = PointLevel(prev_points, prev_diams)

                next_points = np.insert(section.points[dists > distance], 0, bif_pt, axis=0)
                next_diams = np.insert(section.diameters[dists > distance], 0, bif_diam)
                next_pts = PointLevel(next_points, next_diams)

                oblique_pts = PointLevel(oblique_points + next_points[0], oblique_diameters)

                if section.is_root:
                    prev_sec = morphology.append_root_section(prev_pts, APICAL_SUBTYPE_MAP["trunk"])
                else:
                    prev_sec = section.parent.append_section(prev_pts, APICAL_SUBTYPE_MAP["trunk"])
                prev_sec.append_section(next_pts, APICAL_SUBTYPE_MAP["trunk"])
                prev_sec.append_section(oblique_pts, APICAL_SUBTYPE_MAP["oblique"])
                next_sec = prev_sec.children[0]

                _copy(next_sec, section)
                break
    morphology.delete_section(section)
