from pathlib import Path
from mock import patch
from nose.tools import ok_, assert_equal

import pandas as pd
from pandas.testing import assert_frame_equal

import morph_tool.utils as tested

DATA = Path(__file__).resolve().parent / 'data'

def test_is_morphology():
    ok_(tested.is_morphology('a.swc'))
    ok_(tested.is_morphology('a.asc'))
    ok_(tested.is_morphology('a.h5'))
    ok_(not tested.is_morphology('a.blah'))
    ok_(tested.is_morphology('a.blah', extensions={'blah'}))


def test_iter_morphology_files():
    assert_equal(set(tested.iter_morphology_files(DATA / 'folder')),
                 {DATA / 'folder' / 'a.h5',
                  DATA / 'folder' / 'b.swc'})

    assert_equal(set(tested.iter_morphology_files(str(DATA / 'folder'))),
                 {DATA / 'folder' / 'a.h5',
                  DATA / 'folder' / 'b.swc'})

    assert_equal(set(tested.iter_morphology_files(DATA / 'folder', recursive=True)),
                 {DATA / 'folder' / 'a.h5',
                  DATA / 'folder' / 'b.swc',
                  DATA / 'folder' / 'subfolder' / 'g.SWC',
                  DATA / 'folder' / 'subfolder' / 'e.h5',
})

    assert_equal(set(tested.iter_morphology_files(str(DATA / 'folder'), recursive=True)),
                 {DATA / 'folder' / 'a.h5',
                  DATA / 'folder' / 'b.swc',
                  DATA / 'folder' / 'subfolder' / 'g.SWC',
                  DATA / 'folder' / 'subfolder' / 'e.h5',
})

def test_neurondb_dataframe():
    expected = pd.DataFrame(data=[['name1', '1', 'L1_mtype-submtype'],
                                  ['name2', '2', 'L2_bla'],
                                  ['name3', '3', 'L3_tomato']],
                            columns=['name', 'layer', 'mtype'])

    assert_frame_equal(tested.neurondb_dataframe(DATA / 'neurondb.dat'), expected)
