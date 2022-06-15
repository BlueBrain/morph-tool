"""Module to set sub section types, such as oblique for apicals, or collaterals for axons."""
from morphio import SectionType
from neurom import iter_sections
from morphio.mut import Morphology

from morph_tool.apical_point import apical_point_section_segment
from morph_tool.axon_point import axon_point_section


APICAL_SUBTYPE_MAP = {
    "axon": SectionType.axon,
    "basal": SectionType.basal_dendrite,
    "apical": SectionType.apical_dendrite,
    "trunk": SectionType.custom5,
    "oblique": SectionType.custom6,
    "tuft": SectionType.custom7,
}
REV_APICAL_SUBTYPE_MAP = {j: i for i, j in APICAL_SUBTYPE_MAP.items()}

AXON_SUBTYPE_MAP = {
    "axon": SectionType.axon,
    "basal": SectionType.basal_dendrite,
    "apical": SectionType.apical_dendrite,
    "main": SectionType.custom5,
    "collateral": SectionType.custom6,
}
REV_AXON_SUBTYPE_MAP = {j: i for i, j in AXON_SUBTYPE_MAP.items()}


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


def axon_subtype(neuron, axonal_section=None, direction=None, bbox=None, ignore_axis=2):
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


def set_apical_subtypes(neuron, tuft_percent=20):
    """Set apical subtypes to a morphology.

    WARNING: cannot save it to .asc file, use unset_subtype to revert
    """
    subtypes = apical_subtype(neuron, tuft_percent=tuft_percent)
    sections = neuron.sections
    for secid, subtype in subtypes.items():
        sections[secid].morphio_section.type = APICAL_SUBTYPE_MAP[subtype]
    return neuron


def set_axon_subtypes(neuron, **kwargs):
    """Set subtypes to a morphology.

    WARNING: cannot save it to .asc file, use unset_subtype to revert
    """
    subtypes = axon_subtype(neuron, **kwargs)
    sections = neuron.sections
    for secid, subtype in subtypes.items():
        sections[secid].morphio_section.type = AXON_SUBTYPE_MAP[subtype]
    return neuron


def _unset_subtypes(neuron, base_type):
    """Unset subtypes to a morphology.

    Args:
        base_type (SectionType): section type to reset all custom types.

    """
    if isinstance(neuron, Morphology):
        for section in neuron.iter():
            if section.type.name.startswith("custom"):
                section.type = base_type
    else:
        for section in iter_sections(neuron):
            if section.type.name.startswith("custom"):
                section.morphio_section.type = base_type
    return neuron


def unset_apical_subtypes(neuron):
    """Unset apical subtypes to a morphology."""
    return _unset_subtypes(neuron, SectionType.apical_dendrite)


def unset_axon_subtypes(neuron):
    """Unset axon subtypes to a morphology."""
    return _unset_subtypes(neuron, SectionType.axon)
