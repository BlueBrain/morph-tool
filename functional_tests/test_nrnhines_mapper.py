import os
from git import Repo
from tqdm import tqdm


from neurom import load_neuron, COLS
from numpy.testing import assert_almost_equal

from morph_tool.nrnhines import NeuroM_section_to_NRN_section, get_NRN_cell
from morph_tool.utils import iter_morphology_files

from neurom.exceptions import RawDataError

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


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


def _validate_section_mapping(filename, mapping):
    """raise if the mapping NeuroM_section_to_NRN_section is not correct"""
    NeuroM_cell = load_neuron(filename)
    NRN_cell = get_NRN_cell(str(filename))

    for nrm_idx, nrn_idx in mapping.items():
        _validate_section(NRN_cell, NeuroM_cell, nrm_idx, nrn_idx)


def test_mapping_NeuroM_section_to_NRN():
    repo = '/gpfs/bbp.cscs.ch/project/proj30/home/bcoste/MorphologyRepositorySanitized'

    # Concretize list to get tqdm counter
    files = list(iter_morphology_files(repo, recursive=True, extensions={'asc'}))

    for f in tqdm(files):
        try:
            mapping = NeuroM_section_to_NRN_section(f)
            _validate_section_mapping(f, mapping)
        except (RawDataError,
                RuntimeError  # raised by neuron if 2 somas
        ):
            pass
        except:
            print('Error for file: {}'.format(f))
            raise
