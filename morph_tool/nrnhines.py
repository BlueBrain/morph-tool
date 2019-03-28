'''Utils related to the NRN simulator'''
import bluepyopt.ephys as ephys
from numpy.testing import assert_almost_equal
from neurom import COLS, iter_sections, load_neuron
from neurom.core import NeuriteIter
from nose.tools import ok_


def get_cell(filename):
    '''Returns a NRN cell'''
    m = ephys.morphologies.NrnFileMorphology(filename)
    sim = ephys.simulators.NrnSimulator()
    cell = ephys.models.CellModel('test', morph=m, mechs=[])
    cell.instantiate(sim=sim)
    return cell


def _has_single_child(section):
    '''Does the section have only one child ?'''
    return len(section.children()) == 1


def _zero_length_section(section, epsilon=1e-8):
    '''Is zero length section ?'''
    return section.length < epsilon


def _validate_section(nrn_neuron, nrm_neuron, nrm_idx, nrn_idx):
    '''raise if the mapping section_NeuroM_to_NRN is not correct'''
    NRN_sections = list(nrn_neuron.icell.all)

    nrm_s = nrm_neuron.sections[nrm_idx]

    if nrn_idx is None:
        # The only reason for the mapping not to exist is a zero length section
        # because NRN discards them
        ok_(_zero_length_section(nrm_s))
        return

    nrn_s = NRN_sections[nrn_idx]

    nrn_pts = [nrn_s.x3d(0), nrn_s.y3d(0), nrn_s.z3d(0)]
    nrm_pts = nrm_s.points[0, COLS.XYZ]

    err_msg = ('ERROR Section mismatch: NRN ID ({}) != NeuroM ID ({})'
               '{} != {}'.format(nrn_idx, nrm_idx, nrn_pts, nrm_pts))
    assert_almost_equal(nrn_pts,
                        nrm_pts,
                        decimal=2,
                        err_msg=err_msg)


def _validate(filename, mapping):
    '''raise if the mapping section_NeuroM_to_NRN is not correct'''
    cell = get_cell(filename)
    neuron = load_neuron(filename)

    for nrm_idx, nrn_idx in mapping.items():
        _validate_section(cell, neuron, nrm_idx, nrn_idx)


def section_NeuroM_to_NRN(filename):
    '''Returns a mapping from NeuroM section IDs to NRN ones'''
    neurom_cell = load_neuron(filename)

    mapping = dict()

    NRN_sections = list(get_cell(filename).icell.all)

    # Skip soma if exists
    counter = 1 if NRN_sections[0].name().endswith('.soma[0]') else 0

    for NeuroM_section in iter_sections(neurom_cell, neurite_order=NeuriteIter.NRN):
        if _zero_length_section(NeuroM_section):
            mapping[NeuroM_section.id] = None

            if not NeuroM_section.children:
                continue

            NRN_section = NRN_sections[counter]
            counter -= 1

        else:
            mapping[NeuroM_section.id] = counter
            NRN_section = NRN_sections[counter]

        # Skip single child NeuroM_section because they have already been
        # merged in the NeuroM morphology
        while _has_single_child(NRN_section):
            counter += 1
            NRN_section = NRN_section.children()[0]

        counter += 1

    _validate(filename, mapping)
    return mapping
