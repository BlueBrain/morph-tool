
from pathlib import Path
from tempfile import TemporaryDirectory

import morph_tool.morphdb as tested
import numpy as np
import pandas as pd
from nose.tools import assert_raises, eq_
from numpy.testing import assert_array_equal, assert_equal
from pandas.testing import assert_frame_equal

DATA_DIR = Path(__file__).parent / 'data'



def test_from_folder():
    actual = tested.MorphDB.from_folder(DATA_DIR / 'morphdb/from_folder',
                                        mtypes={'simple': 'typeA', 'simple2': 'typeB:withsubtype'},
                                        label='release-1').df

    expected = pd.read_csv(DATA_DIR / 'morphdb/from_folder/expected.csv', header=0, keep_default_na=False, na_values={'path': ['']})
    expected.path = expected.path.astype(str).replace(to_replace={'nan': None})

    assert_frame_equal(actual.drop(columns='axon_inputs'),
                       expected.drop(columns='axon_inputs'),)


def test_from_neurondb():
    actual = tested.MorphDB.from_neurondb(DATA_DIR / 'morphdb/from_neurondb/neurondb.xml',
                                          label='release-1').df

    expected = pd.read_csv(DATA_DIR / 'morphdb/from_neurondb/expected.csv', header=0, keep_default_na=False, na_values={'path': ['']})

    expected.layer = expected.layer.astype(str)
    expected.path = expected.path.astype(str).replace(to_replace={'nan': None})

    assert_frame_equal(actual.drop(columns='axon_inputs'),
                       expected.drop(columns='axon_inputs'),)

def test_single_axon_input():
    actual = tested.MorphDB.from_neurondb(
        DATA_DIR / 'morphdb/from_neurondb/single-axon-input.neurondb',
        label='release-1').df


    assert_equal(actual.axon_inputs.iloc[0], ['C270106A'])


def test_MorphInfo():
    morph = tested.MorphInfo(name='a', mtype='b', layer='c')
    assert_equal(str(morph), 'MorphInfo(name=a, mtype=b, layer=c)')


def test_read_msubtype():
    morphology_folder = DATA_DIR / 'morphdb/from_neurondb/'
    df = tested.MorphDB.from_neurondb(morphology_folder / 'neurondb-msubtype.xml').df
    columns = ['mtype', 'msubtype', 'mtype_no_subtype']
    assert_frame_equal(df[columns], pd.DataFrame(data=[['L1_DAC:A', 'A', 'L1_DAC'],
                                                       ['L1_DAC', '', 'L1_DAC'],
                                                       ['L1_DAC', '', 'L1_DAC']],
                                                 columns=columns))



def test_write_neurondb_dat():
    morphology_folder = DATA_DIR / 'morphdb/from_neurondb/'
    original = tested.MorphDB.from_neurondb(morphology_folder / 'neurondb-msubtype.xml')

    with TemporaryDirectory() as temp_dir:
        path = Path(temp_dir, f'neurondb.dat')
        original.write(path)

        new = tested.MorphDB.from_neurondb(path, morphology_folder=morphology_folder)
        assert_frame_equal(original.df, new.df)


def test_add():
    original = tested.MorphDB.from_neurondb(DATA_DIR / 'morphdb/from_neurondb/neurondb-msubtype.xml')
    morphs = [tested.MorphInfo(name='tomato', mtype='banana:split', layer=1),
              tested.MorphInfo(name='elon', mtype='musk', layer='2')]

    # testing adding lists
    total = original + morphs
    assert_array_equal(total.df.name, ['tkb061126a4_ch0_cc2_h_zk_60x_1', 'missing-morph', 'simple3',
       'tomato', 'elon'])

    # testing adding an MorphDB
    total = original + original
    assert_array_equal(total.df.name, ['tkb061126a4_ch0_cc2_h_zk_60x_1', 'missing-morph', 'simple3',
       'tkb061126a4_ch0_cc2_h_zk_60x_1', 'missing-morph', 'simple3'])

    assert_raises(TypeError, original.__add__, None)



def test_iadd():
    original = tested.MorphDB.from_neurondb(DATA_DIR / 'morphdb/from_neurondb/neurondb-msubtype.xml')
    morphs = [tested.MorphInfo(name='tomato', mtype='banana:split', layer=1),
              tested.MorphInfo(name='elon', mtype='musk', layer='2')]

    original += morphs
    assert_array_equal(original.df.name,
                       ['tkb061126a4_ch0_cc2_h_zk_60x_1', 'missing-morph', 'simple3',
                        'tomato', 'elon'])


    original = tested.MorphDB.from_neurondb(DATA_DIR / 'morphdb/from_neurondb/neurondb-msubtype.xml')
    original += original
    assert_array_equal(original.df.name, ['tkb061126a4_ch0_cc2_h_zk_60x_1', 'missing-morph', 'simple3',
       'tkb061126a4_ch0_cc2_h_zk_60x_1', 'missing-morph', 'simple3'])

    assert_raises(TypeError, original.__iadd__, None)

def test_write_neurondb_xml():
    morphology_folder = DATA_DIR / 'morphdb/from_neurondb/'
    original = tested.MorphDB.from_neurondb(morphology_folder / 'neurondb-msubtype.xml')

    with TemporaryDirectory() as temp_dir:
        path = Path(temp_dir, 'neurondb.xml')
        original.write(path)

        new = tested.MorphDB.from_neurondb(path, morphology_folder=morphology_folder)
        assert_frame_equal(original.df, new.df)

def test_load_raises():
    original = tested.MorphDB.from_neurondb(
        DATA_DIR / 'morphdb/from_neurondb/neurondb-only-dat-info.xml')
    assert_raises(ValueError, original.write, Path('neurondb.wrong-format'))


def test_features():
    original = tested.MorphDB.from_neurondb(
        DATA_DIR / 'morphdb/from_neurondb/neurondb-only-dat-info.xml')
    assert_raises(ValueError, original.features,
                  {'neurite': {'section_lengths': ['max']}})

    original.df = original.df[~original.df.path.isnull()]
    features = original.features({'neurite': {'section_lengths': ['max']}})
    # features.to_csv(DATA_DIR / 'morphdb/from_neurondb/features.csv', index=False)
    expected = pd.read_csv(DATA_DIR / 'morphdb/from_neurondb/features.csv',
                           header=[0, 1])
    expected['neuron'] = expected['neuron'].fillna('')
    expected['neuron', 'layer'] = expected['neuron', 'layer'].astype(str)
    expected['neuron', 'path'] = expected['neuron', 'path'].apply(Path)
    for key in tested.BOOLEAN_REPAIR_ATTRS:
        expected['neuron', key] = expected['neuron', key].astype(bool)
    assert_frame_equal(features.drop(columns=('neuron', 'axon_inputs')),
                       expected.drop(columns=('neuron', 'axon_inputs')))
