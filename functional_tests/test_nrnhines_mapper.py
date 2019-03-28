import os
from git import Repo
from tqdm import tqdm


from morph_tool.nrnhines import section_NeuroM_to_NRN
from morph_tool.utils import iter_morphology_files

from neurom.exceptions import RawDataError

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

def test_mapping_section_NeuroM_to_NRN():
    repo = os.path.join('/tmp', 'MorphologyRepository')
    if not os.path.exists(repo):
        Repo.clone_from('ssh://bbpcode.epfl.ch/experiment/MorphologyRepository',
                        repo, depth=1)

    repo = '/gpfs/bbp.cscs.ch/project/proj68/tmp/NCX-83/one-column-proj64/20190308/synthesis/morphologies/hashed/06'

    for f in tqdm(iter_morphology_files(repo, extensions={'asc'})):
        try:
            section_NeuroM_to_NRN(f)
        except RawDataError:
            pass
