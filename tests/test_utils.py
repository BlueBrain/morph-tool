from mock import patch
from nose.tools import ok_, assert_equal

from morph_tool.utils import is_morphology, iter_morphology_files

def test_is_morphology():
    ok_(is_morphology('a.swc'))
    ok_(is_morphology('a.asc'))
    ok_(is_morphology('a.h5'))
    ok_(not is_morphology('a.blah'))
    ok_(is_morphology('a.blah', extensions={'blah'}))


@patch('morph_tool.utils.os.listdir', return_value=['a.h5', 'b.swc', 'c.blah', 'd'])
@patch('morph_tool.utils.os.walk', return_value=[('folder', 'subfolders', ['a.h5', 'c.blah']),
                                                 ('folder/subfolder', 'subfolders', ['b.swc'])])
def test_iter_morphology_files(_, __):
    assert_equal(list(iter_morphology_files('folder')), ['folder/a.h5', 'folder/b.swc'])
    assert_equal(list(iter_morphology_files('folder', recursive=True)),
                 ['folder/a.h5', 'folder/subfolder/b.swc'])
