from pathlib2 import Path
from mock import patch
from nose.tools import ok_, assert_equal

from morph_tool.utils import is_morphology, iter_morphology_files

PATH = Path(__file__).resolve().parent / 'data'

def test_is_morphology():
    ok_(is_morphology('a.swc'))
    ok_(is_morphology('a.asc'))
    ok_(is_morphology('a.h5'))
    ok_(not is_morphology('a.blah'))
    ok_(is_morphology('a.blah', extensions={'blah'}))


def test_iter_morphology_files():
    assert_equal(set(iter_morphology_files(PATH / 'folder')),
                 {PATH / 'folder' / 'a.h5',
                  PATH / 'folder' / 'b.swc'})

    assert_equal(set(iter_morphology_files(str(PATH / 'folder'))),
                 {PATH / 'folder' / 'a.h5',
                  PATH / 'folder' / 'b.swc'})

    assert_equal(set(iter_morphology_files(PATH / 'folder', recursive=True)),
                 {PATH / 'folder' / 'a.h5',
                  PATH / 'folder' / 'b.swc',
                  PATH / 'folder' / 'subfolder' / 'g.SWC',
                  PATH / 'folder' / 'subfolder' / 'e.h5',
})

    assert_equal(set(iter_morphology_files(str(PATH / 'folder'), recursive=True)),
                 {PATH / 'folder' / 'a.h5',
                  PATH / 'folder' / 'b.swc',
                  PATH / 'folder' / 'subfolder' / 'g.SWC',
                  PATH / 'folder' / 'subfolder' / 'e.h5',
})
