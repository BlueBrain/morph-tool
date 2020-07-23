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

    expected_axon_input = [
        'C270106A',
        'C270106C',
        'C270106G',
        'sm080529a1-5_idA',
        'sm080619a1-7_idG',
        'sm080625a1-6_idD',
        'sm080723a1-4_idC',
        'sm080930a1-5_idC',
    ]

    expected = pd.DataFrame([row for row
                             in [['C270106A', 'L1_DAC', '', 'L1_DAC', '1', True, True, False, False, True, True, True, True, True, expected_axon_input],
                                 ['C270106C', 'L1_DAC', '', 'L1_DAC', '1', True, True, True, True, True, True, True, True, True, expected_axon_input],
                                 ['a_neuron', 'an_mtype', 'a_subtype', 'an_mtype:a_subtype', '1', False, True, True, True, True, True, True, True, True, []],
                                 ['a_2nd_neuron', 'an_mtype', 'a_subtype', 'an_mtype:a_subtype', '1', True, True, True, True, True, True, True, True, True, []],
                                 ]],
                            columns=['name', 'mtype', 'msubtype', 'fullmtype', 'layer',
                                     'use_axon', 'use_dendrites', 'axon_repair', 'dendrite_repair', 'basal_dendrite_repair',
                                     'tuft_dendrite_repair', 'oblique_dendrite_repair', 'unravel',
                                     'use_for_stats', 'axon_inputs'])

    assert_frame_equal(df, expected)

    assert_raises(ValueError, tested.neurondb_dataframe, DATA / 'neurondb.wrongext')
