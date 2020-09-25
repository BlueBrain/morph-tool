import sys
from io import StringIO
from pathlib import Path

import pandas as pd
from mock import patch
from nose.tools import assert_equal, assert_raises, ok_
from numpy.testing import assert_array_equal
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

def test_find_morph():
    folder = DATA / 'test-neurondb-with-path'
    assert_equal(tested.find_morph(folder, 'not-here.h5'),
                 None)
    assert_equal(tested.find_morph(folder, 'C270106A'),
                 folder / 'C270106A.h5')
    assert_equal(tested.find_morph(folder, 'C270106C.wrongext'),
                 None)


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


def test_neurondb_dataframe_no_repair():
    df = tested.neurondb_dataframe(DATA / 'neurondb-no-repair.xml')
    expected = pd.DataFrame(data=[['C270106A', '1', 'L1_DAC', True],
                                  ['C270106C', '1', 'L1_DAC', True],

                                  # This time (cf test_neurondb_dataframe), use_axon is True
                                  ['a_neuron', '1', 'an_mtype:a_subtype', True],

                                  ['a_2nd_neuron', '1', 'an_mtype:a_subtype', True],
                                  ],
                            columns=['name', 'layer', 'mtype', 'use_axon'])

    assert_frame_equal(df, expected)


def mock_path_content(content):
    class MockPathContent:
        '''pathlib.Path mock to mock the call to Path.open()'''
        def __enter__(self):
            return StringIO(content)

        def __exit__(self, exc_type, exc_val, exc_tb):
            pass
    return MockPathContent


def test_neurondb_dataframe_single_morph():
    neurondb = DATA / 'neurondb.xml'

    neurondb_template = '''
    <neurondb>
    <listing>
    <morphology>
      <name>C270106A</name>
      <mtype>L1_DAC</mtype>
      <msubtype></msubtype>
      <layer>1</layer>
      <repair>
        <use_axon>True</use_axon>
      </repair>
    </morphology>
    </listing>
    </neurondb>
    '''

    with patch.object(Path, 'open', mock_path_content(neurondb_template)):
        df = tested.neurondb_dataframe(neurondb)
        expected = pd.DataFrame(data=[['C270106A', '1', 'L1_DAC', True]],
                                columns=['name', 'layer', 'mtype', 'use_axon'])
        assert_frame_equal(df, expected)

def test_neurondb_dataframe_use_axon():
    neurondb = DATA / 'neurondb.xml'

    neurondb_template = '''
    <neurondb>
    <listing>
    <morphology>
      <name>C270106A</name>
      <mtype>L1_DAC</mtype>
      <msubtype></msubtype>
      <layer>1</layer>
      <repair>
        <use_axon>{}</use_axon>
      </repair>
    </morphology>
    </listing>
    </neurondb>
    '''

    for use_axon in ['True', 'true', '']:
        with patch.object(Path, 'open', mock_path_content(neurondb_template.format(use_axon))):
            df = tested.neurondb_dataframe(neurondb)
            assert_equal(df.loc[0, 'use_axon'], True)

    for use_axon in ['False', 'false']:
        with patch.object(Path, 'open', mock_path_content(neurondb_template.format(use_axon))):
            df = tested.neurondb_dataframe(neurondb)
            assert_equal(df.loc[0, 'use_axon'], False)

    for use_axon in ['tRuE', 'fals', 0, 1, 'mickael jackson']:
        with patch.object(Path, 'open', mock_path_content(neurondb_template.format(use_axon))):
            assert_raises(AssertionError, tested.neurondb_dataframe, neurondb)


def test_neurondb_dataframe_with_path():
    folder = DATA / 'test-neurondb-with-path'
    df = tested.neurondb_dataframe(DATA / 'neurondb.xml', morphology_dir=folder)

    assert_array_equal([p.stem if p is not None else None for p in df.path],
                       ['C270106A', None, 'a_neuron', None])
