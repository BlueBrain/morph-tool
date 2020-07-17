from pathlib import Path

import pandas as pd
from nose.tools import assert_equal, assert_raises, ok_
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
    df = tested.neurondb_dataframe(DATA / 'neurondb.xml')
    expected = pd.DataFrame(data=[['C270106A', '1', 'L1_DAC', True],
                                  ['C270106C', '1', 'L1_DAC', True],
                                  ['a_neuron', '1', 'an_mtype:a_subtype', False],
                                  ['a_2nd_neuron', '1', 'an_mtype:a_subtype', True],
                                  ],
                            columns=['name', 'layer', 'mtype', 'use_axon'])

    assert_frame_equal(df, expected)

    assert_raises(ValueError, tested.neurondb_dataframe, DATA / 'neurondb.wrongext')
