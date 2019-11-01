"""Utils related to the NRN simulator"""
import logging
import numpy as np
from numpy.testing import assert_almost_equal
from neurom import COLS, iter_sections, load_neuron, NeuriteType
from neurom.core import NeuriteIter

L = logging.getLogger('morph_tool')


def _get_NRN_cell(filename):
    """Returns a NRN cell"""
    try:
        # pylint: disable=import-outside-toplevel
        import bluepyopt.ephys as ephys
    except ImportError:
        raise Exception('_get_NRN_cell requires the extra: all. '
                        'Please install it by doing: pip install morph-tool[all]')

    m = ephys.morphologies.NrnFileMorphology(filename)
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


def _validate_section(nrn_neuron, nrm_neuron, nrm_idx, nrn_idx):
    """raise if the mapping NeuroM_section_to_NRN_section is not correct"""
    NRN_sections = list(nrn_neuron.icell.all)

    nrm_s = nrm_neuron.sections[nrm_idx]

    if nrn_idx is None:
        # The only reason for the mapping not to exist is a zero length section
        # because NRN discards them
        assert _zero_length_section(nrm_s)
        return

    nrn_s = NRN_sections[nrn_idx]

    nrn_pts = [nrn_s.x3d(0), nrn_s.y3d(0), nrn_s.z3d(0)]
    nrm_pts = nrm_s.points[0, COLS.XYZ]

    err_msg = ('ERROR Section mismatch: NRN ID ({}) != NeuroM ID ({})'
               'NRN section:\n{}\n\nNeuroM section:\n{}'.format(
                   nrn_idx, nrm_idx,
                   nrn_pts, nrm_pts))
    assert_almost_equal(nrn_pts,
                        nrm_pts,
                        decimal=2,
                        err_msg=err_msg)


def _validate_section_mapping(NeuroM_cell, NRN_cell, mapping):
    """raise if the mapping NeuroM_section_to_NRN_section is not correct"""
    for nrm_idx, nrn_idx in mapping.items():
        _validate_section(NRN_cell, NeuroM_cell, nrm_idx, nrn_idx)


def NeuroM_section_to_NRN_section(filename):
    """Returns a mapping from NeuroM section IDs to NRN ones"""
    NeuroM_cell = load_neuron(filename)
    NRN_cell = _get_NRN_cell(filename)

    mapping = dict()

    NRN_sections = list(NRN_cell.icell.all)

    def is_soma(NRN_section):
        """Is the NRN section a soma section"""
        return NRN_section.name().endswith('.soma[0]')

    # Skip soma if exists
    counter = 1 if is_soma(NRN_sections[0]) else 0

    for NeuroM_section in iter_sections(NeuroM_cell, neurite_order=NeuriteIter.NRN):
        if _zero_length_section(NeuroM_section):
            mapping[NeuroM_section.id] = None

            if not NeuroM_section.children:
                L.debug('Zero length section without children (NeuroM section id: %s)',
                        NeuroM_section.id)
                continue

            L.debug('Zero length section with children')
            NRN_section = NRN_sections[counter]
            counter -= 1

        else:
            mapping[NeuroM_section.id] = counter
            NRN_section = NRN_sections[counter]

        L.debug('NeuroM section (%s) has been mapped to NRN section (%s)',
                NeuroM_section.id, mapping[NeuroM_section.id])

        # Skip single child NeuroM_section because they have already been
        # merged in the NeuroM morphology
        while _has_single_child(NRN_section):
            L.debug('Skipping single child')
            counter += 1
            NRN_section = NRN_section.children()[0]

        counter += 1

    _validate_section_mapping(NeuroM_cell, NRN_cell, mapping)
    return mapping


def _interpolate_compartments(points, boundaries_segment_ids, boundaries_positions):
    """
    Returns the path of each compartment based on its starting and ending position

    Args:
        points (numpy.array): a 3D array of the section points
        boundaries_segment_ids (List): the ids of the segment in which each
                                       compartment boundary belongs
        boundaries_positions (List): the 3D positions of the compartment boundaries
    """
    compartment_points = list()
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
    """
    Returns the list of paths (list of 3D points) that form each compartment
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
        cumulative_pathlength, pathlengths_at_compartment_boundaries, side='right') - 1).tolist()
    boundaries_positions = list()
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


def NeuroM_section_to_NRN_compartment_paths(morph_path):
    """Returns a dictionary NeuroM section id -> path of each compartment for the section

    Path are formed by following the section points until the pathlength of the compartment is
    reached.

    Args:
        morph_path (str): the morphology path


    1) Compute the cumulative pathlength along the section segments_directions
    2) Get the compartment pathlengths (compartment are of equal pathlength in a given section)
    3) Compartment by compartment, follow the points until the compartment pathlength is reached


    Example for one section:

                   (1, 2) ------ (2, 2)
                      |
                      |
                      |
    (0, 0) ------- (1, 0)


    If n_compartments == 3, three paths are returned:
    [array([[0.        , 0.        , 0.        ],
            [1.        , 0.        , 0.        ],
            [1.        , 0.33333333, 0.        ]]),

     array([[1.        , 0.33333333, 0.        ],
            [1.        , 1.66666667, 0.        ]]),

     array([[1.        , 1.66666667, 0.        ],
            [1.        , 2.        , 0.        ],
            [2.        , 2.        , 0.        ]])]
    """

    NeuroM_cell = load_neuron(morph_path)
    NRN_neuron = _get_NRN_cell(morph_path)
    NRN_sections = list(NRN_neuron.icell.all)

    mapping = NeuroM_section_to_NRN_section(morph_path)

    NeuroM_to_compartment_position_mapping = dict()

    for section in NeuroM_cell.sections:
        if section.type == NeuriteType.soma:
            continue

        NRN_section = NRN_sections[mapping[section.id]]

        NeuroM_to_compartment_position_mapping[section.id] = _compartment_paths(
            section.points[:, COLS.XYZ], NRN_section.nseg)

    return NeuroM_to_compartment_position_mapping
