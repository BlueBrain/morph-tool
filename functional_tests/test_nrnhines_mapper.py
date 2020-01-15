import os
from git import Repo
from tqdm import tqdm


from morph_tool.nrnhines import NeuroM_section_to_NRN_section
from morph_tool.utils import iter_morphology_files

from neurom.exceptions import RawDataError

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

def test_mapping_NeuroM_section_to_NRN():
    repo = os.path.join('/tmp', 'MorphologyRepository')
    if not os.path.exists(repo):
        Repo.clone_from('ssh://bbpcode.epfl.ch/experiment/MorphologyRepository',
                        repo, depth=1)

    repo = '/gpfs/bbp.cscs.ch/project/proj68/tmp/NCX-83/one-column-proj64/20190308/synthesis/morphologies/hashed/06'

    # Concretize list to get tqdm counter
    files = list(iter_morphology_files(repo, recursive=True, extensions={'asc'}))

    for f in tqdm(files):
        try:
            NeuroM_section_to_NRN_section(str(f))
        except (RawDataError,
                RuntimeError  # raised by neuron if 2 somas
        ):
            pass
        except:
            print('Error for file: {}'.format(f))
            raise
