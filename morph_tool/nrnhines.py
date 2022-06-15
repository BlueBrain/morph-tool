"""Utils related to the NRN simulator."""
import logging
import multiprocessing.pool
from pathlib import Path
from typing import List, Sequence, Union

import numpy as np
from neurom import COLS, NeuriteType, iter_sections, load_morphology
from neurom.core.types import NeuriteIter
from numpy.testing import assert_almost_equal

try:
    import neuron
except ImportError as e:
    raise ImportError(
        'morph-tool[nrn] is not installed. Please install: pip install morph-tool[nrn]'
    ) from e

L = logging.getLogger(__name__)


def get_NRN_cell(filename: Path):
    """Returns a NRN cell."""
    try:
        # pylint: disable=import-outside-toplevel
        from bluepyopt import ephys
    except ImportError as e_:
        raise ImportError(
            'bluepyopt not installed; please use `pip install morph-tool[nrn]`') from e_
    m = ephys.morphologies.NrnFileMorphology(str(filename))
    sim = ephys.simulators.NrnSimulator()
    cell = ephys.models.CellModel('test', morph=m, mechs=[])
    cell.instantiate(sim=sim)
    return cell


def _has_single_child(section):
    """Does the section have only one child ?"""
    return len(section.children()) == 1


def _zero_length_section(section, epsilon=1e-8):
    """Is zero length section ?"""
    return section.length < epsilon


def _validate_section(nrn_neuron, nrm_neuron, nrm_idx, nrn_idx, neurite_type):
    """Raise if the mapping NeuroM_section_to_NRN_section is not correct."""
    nrn_sections = list(getattr(nrn_neuron.icell, neurite_type.split("_")[0]))

    if neurite_type != "all":
        nrn_section_counter = 0
        for nrm_s in iter_sections(nrm_neuron, neurite_order=NeuriteIter.NRN):
            if _has_type(nrm_s, neurite_type):
                if nrn_section_counter == nrn_idx:
                    break
                nrn_section_counter += 1
    else:
        nrm_s = nrm_neuron.sections[nrm_idx]

    if nrn_idx is None:
        # The only reason for the mapping not to exist is a zero length section
        # because NRN discards them
        assert _zero_length_section(nrm_s)
        return

    nrn_s = nrn_sections[nrn_idx]

    nrn_pts = [nrn_s.x3d(0), nrn_s.y3d(0), nrn_s.z3d(0)]
    nrm_pts = nrm_s.points[0, COLS.XYZ]

    err_msg = (f'ERROR Section mismatch: NRN ID ({nrn_idx}) != NeuroM ID ({nrm_idx})'
               'NRN section:\n'
               f'{nrn_pts}\n\n'
               'NeuroM section:\n'
               f'{nrm_pts}')
    assert_almost_equal(nrn_pts,
                        nrm_pts,
                        decimal=2,
                        err_msg=err_msg)


def _validate_section_mapping(neurom_cell, nrn_cell, mapping, neurite_type):
    """Raise if the mapping NeuroM_section_to_NRN_section is not correct."""
    for nrm_idx, nrn_idx in mapping.items():
        _validate_section(nrn_cell, neurom_cell, nrm_idx, nrn_idx, neurite_type)


def _has_type(section, neurite_type):
    """Check if the NeuriteType of section is the expected one."""
    if neurite_type == "all":
        return True
    return section.type == getattr(NeuriteType, neurite_type)


def NeuroM_section_to_NRN_section(filename: Path, neurite_type: str = "all", reverse: bool = False):
    """Returns a mapping from NeuroM section IDs to NRN ones.

    In NEURON, neurites are implicitely defined by lists of sections,
    and the index of the sections in each list starts from 0. Thus the index
    of a section  in the list of apical sections will differ from the same section
    in the list of all sections. With the argument `neurite_type` one can match ids
    of sections whithin each NEURON list of sections.

    Args:
        filename: path to morphology file
        neurite_type: type of neurite to convert section ids (should be NeuriteType name)
        reverse: returns the reverse mapping, between NRN to NeuroM
    """
    assert hasattr(NeuriteType, neurite_type), f"{neurite_type} is not a valid NeuriteType"

    neurom_cell = load_morphology(filename)
    nrn_cell = get_NRN_cell(filename)
    nrn_sections = list(getattr(nrn_cell.icell, neurite_type.split("_")[0]))

    def is_soma(nrn_section):
        """Is the NRN section a soma section."""
        return nrn_section.name().endswith(".soma[0]")

    # Skip soma if exists
    nrn_section_counter = 1 if is_soma(nrn_sections[0]) else 0

    nrm_to_nrn_sections = {}
    for neurom_section in iter_sections(neurom_cell, neurite_order=NeuriteIter.NRN):
        if _has_type(neurom_section, neurite_type):
            if _zero_length_section(neurom_section):
                nrm_to_nrn_sections[neurom_section.id] = None

                if not neurom_section.children:
                    L.debug(
                        "Zero length section without children (NeuroM section id: %s)",
                        neurom_section.id,
                    )
                    continue

                L.debug("Zero length section with children")
                nrn_section_counter -= 1

            else:
                nrm_to_nrn_sections[neurom_section.id] = nrn_section_counter

            L.debug('NeuroM section (%s) has been mapped to NRN section (%s)',
                    neurom_section.id, nrm_to_nrn_sections[neurom_section.id])
            nrn_section_counter += 1

    _validate_section_mapping(neurom_cell, nrn_cell, nrm_to_nrn_sections, neurite_type)
    return {j: i for i, j in nrm_to_nrn_sections.items()} if reverse else nrm_to_nrn_sections


def _interpolate_compartments(points, boundaries_segment_ids, boundaries_positions):
    """Returns the path of each compartment based on its starting and ending position.

    Args:
        points (numpy.array): a 3D array of the section points
        boundaries_segment_ids (List): the ids of the segment in which each
                                       compartment boundary belongs
        boundaries_positions (List): the 3D positions of the compartment boundaries
    """
    compartment_points = []
    for i, (segment_id_start, position) in enumerate(
            zip(boundaries_segment_ids[:-1], boundaries_positions[:-1])):
        compartment = [position]
        segment_id_end = boundaries_segment_ids[i + 1]

        # The compartment might span on multiple segments, we add the intermediate points here
        intermediate_point_id = segment_id_start + 1
        while intermediate_point_id <= segment_id_end:
            compartment.append(points[intermediate_point_id])
            intermediate_point_id += 1

        # If the boundary point matches a point from the 'points' list, it has already been added
        # to the 'compartment' list, no need to re-add it
        if not np.allclose(compartment[-1], boundaries_positions[i + 1]):
            compartment.append(boundaries_positions[i + 1])

        compartment_points.append(np.vstack(compartment))
    return compartment_points


def _compartment_paths(points, n_compartments):
    """Returns the list of paths (list of 3D points) that form each compartment.

    Args:
        points (numpy.array): a 3D array of points
        n_compartments (int): the number of compartments
    """
    segments_directions = np.diff(points, axis=0)
    segment_lengths = np.linalg.norm(segments_directions, axis=1)
    cumulative_pathlength = np.append(0, np.cumsum(segment_lengths))
    pathlengths_at_compartment_boundaries = [
        i * cumulative_pathlength[-1] / float(n_compartments) for i in range(n_compartments)]

    boundaries_segment_ids = (np.searchsorted(
        cumulative_pathlength, pathlengths_at_compartment_boundaries, side="right") - 1).tolist()
    boundaries_positions = []
    for segment_id, boundary_pathlength in zip(boundaries_segment_ids,
                                               pathlengths_at_compartment_boundaries):
        # The pathlength between the start of segment #segment_id and the boundary
        remaining_pathlength = boundary_pathlength - cumulative_pathlength[segment_id]

        # the boundary is somewhere in between point #segment_id and #segment_id+1
        # Here we compute its position in term of relative pathlength
        segment_fraction = (remaining_pathlength / segment_lengths[segment_id])
        position = points[segment_id] + segment_fraction * segments_directions[segment_id]
        boundaries_positions.append(position)

    # Adding the last boundary which corresponds to the last point of the section
    boundaries_segment_ids.append(len(points) - 1)
    boundaries_positions.append(points[-1])

    return _interpolate_compartments(points, boundaries_segment_ids, boundaries_positions)


def NeuroM_section_to_NRN_compartment_paths(morph_path: Path):
    """Returns a dictionary NeuroM section id -> path of each compartment for the section.

    Path are formed by following the section points until the pathlength of the compartment is
    reached.

    Args:
        morph_path: the morphology path

    1) Compute the cumulative pathlength along the section segments_directions
    2) Get the compartment pathlengths (compartment are of equal pathlength in a given section)
    3) Compartment by compartment, follow the points until the compartment pathlength is reached


    Example for one section::

                       (1, 2) ------ (2, 2)
                          |
                          |
                          |
        (0, 0) ------- (1, 0)


    If n_compartments == 3, three paths are returned::

        [array([[0.        , 0.        , 0.        ],
                [1.        , 0.        , 0.        ],
                [1.        , 0.33333333, 0.        ]]),

         array([[1.        , 0.33333333, 0.        ],
                [1.        , 1.66666667, 0.        ]]),

         array([[1.        , 1.66666667, 0.        ],
                [1.        , 2.        , 0.        ],
                [2.        , 2.        , 0.        ]])]
    """
    neurom_cell = load_morphology(morph_path)
    nrn_neuron = get_NRN_cell(morph_path)
    nrn_sections = list(nrn_neuron.icell.all)

    mapping = NeuroM_section_to_NRN_section(morph_path)

    neurom_to_compartment_position_mapping = {}
    for section in neurom_cell.sections:
        if section.type == NeuriteType.soma:
            continue

        nrn_section = nrn_sections[mapping[section.id]]

        neurom_to_compartment_position_mapping[section.id] = _compartment_paths(
            section.points[:, COLS.XYZ], nrn_section.nseg
        )

    return neurom_to_compartment_position_mapping


def point_to_section_end(sections: Sequence[neuron.nrn.Section],  # pylint: disable=no-member
                         point: List[float],
                         atol: float = 1e-08,
                         rtol: float = 1e-05) -> Union[None, int]:
    """Returns the index of the first section whose end is close to ``point``.

    The distance between the end and ``point`` must be less than EPSILON from POINT. If no section
    satisfies this requirement, returns None. The iteration order is given by the section index.

    Args:
        sections: a sequence of sections
        point: the points's 3D coordinates
        atol: absolute tolerance
        rtol: relative tolerance
    """
    point = np.asarray(point)

    for index, section in enumerate(sections):
        last_index = section.n3d() - 1
        last_section_point = [section.x3d(last_index),
                              section.y3d(last_index),
                              section.z3d(last_index)]

        if np.isclose(point, last_section_point, atol=atol, rtol=rtol).all():
            return index
    return None


class NestedPool(multiprocessing.pool.Pool):  # pylint: disable=abstract-method
    """Class that represents a MultiProcessing nested pool."""

    class Process(multiprocessing.Process):
        """Class that represents a non-daemon process."""

        def __init__(self, *args, **kwargs):
            """Ensures group=None to avoid exception from multiprocessing.

            Note:
                The exception raised by multiprocessing is:
                AssertionError: group argument must be None for now
            """
            kwargs.pop("group", None)
            super().__init__(None, *args[1:], **kwargs)

        @property
        def daemon(self):
            """Stub method. Always returns False."""
            return False

        @daemon.setter
        def daemon(self, value):
            """Stub method. Does not do anything."""


def isolate(func):
    """Isolate a generic function for independent NEURON instances.

    It must be used in conjunction with NestedPool.

    Example:

    .. code-block:: python

        def _to_be_isolated(morphology_path, point):
            cell = nrnhines.get_NRN_cell(morphology_path)
            return nrnhines.point_to_section_end(cell.icell.all, point)

        def _isolated(morph_data):
            return nrnhines.isolate(_to_be_isolated)(*morph_data)

        with nrnhines.NestedPool(processes=n_workers) as pool:
            result = pool.imap_unordered(_isolated, data)


    Args:
        func (function): function to isolate

    Returns:
        the isolated function

    Note: it does not work as decorator.
    """
    def func_isolated(*args, **kwargs):
        with NestedPool(1, maxtasksperchild=1) as pool:
            return pool.apply(func, args, kwargs)
    return func_isolated
