import csv
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from nose.tools import eq_, ok_, assert_not_equal
from numpy.testing import assert_array_equal

from morph_tool import morphdb, utils
from morph_tool.morphdb import MorphInfo, MorphologyDB

# from .utils import create_fake_files, TemporaryDirectory, write_neurondb

DATA_PATH = Path(__file__).parent / 'data'


def test_MorphInfo_hash():
    eq_(hash(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='1'))),
        hash(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='1', use_axon=True))),)

    assert_not_equal(hash(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='1'))),
                     hash(MorphInfo(dict(name='morph0', mtype='L1_fake_type0:A', layer='1'))),)

    assert_not_equal(hash(MorphInfo(dict(name='aaa', mtype='L1_fake_type0', layer='1'))),
                     hash(MorphInfo(dict(name='bbb', mtype='L1_fake_type0', layer='1'))),)

    assert_not_equal(hash(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='1'))),
                     hash(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='2'))),)


def test__xmlmorphinfo_from_xml():
    expected_list = [['sm080529a1-5_idB', 'L1_SAC', '1', True],
                     ['sm080529a1-5_idB', 'L1_SAC', '1', False],
                     ['sm080529a1-5_idB', 'L1_SAC', '1', True]]
    for morph_info, expected in zip(MorphologyDB(DATA_PATH / 'neurondb2.xml'),
                                    expected_list):
        eq_(morph_info.name, expected[0])
        eq_(morph_info.mtype, expected[1])
        eq_(morph_info.layer, expected[2])
        eq_(morph_info.use_axon, expected[3])


def test_remove_morphs():
    with TemporaryDirectory() as temp_dir:
        morph_db = MorphologyDB(DATA_PATH / 'neurondb3.xml')
        morph_db.add_morph(MorphInfo(dict(name='fake_name', mtype='fake_mtype', layer='0')))
        morph_db.remove_morphs(('fake_name', 'C270106A', ))
        names = [morph.name for morph in morph_db.morphologies]
        ok_('fake_name' not in names)
        ok_('C270106A' not in names)

    # make sure that removing morphs correctly sets the _known_morphology set
    morph_db = morphdb.MorphologyDB()
    morph_db.add_morph(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='1')))
    morph_db.add_morph(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='2')))
    morph_db.remove_morphs(('morph0', ))
    eq_(0, len(list(morph_db)))
    ok_(morph_db.add_morph(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='1'))))
    ok_(morph_db.add_morph(MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer='2'))))


def test_load_neurondb_xml():
    neurondb_xml = Path(DATA_PATH, 'neurondb3.xml')
    morph_db = MorphologyDB(neurondb_xml)
    eq_(len(morph_db.morphologies), 6)


def test_write_neurondb_dat():
    neurondb_xml = Path(DATA_PATH, 'neurondb3.xml')
    morph_db = MorphologyDB(neurondb_xml)
    with TemporaryDirectory() as temp_dir:
        for ext in ['dat', 'csv']:
            path = Path(temp_dir, f'neurondb.{ext}')
            morph_db.write(path)
            with path.open() as fd:
                lines = fd.readlines()
            if ext == 'dat':
                lines = [line.split() for line in lines]
            else:
                # first line is column names, hence is skipped
                lines = [line.strip().split(',') for line in lines[1:]]

            eq_(lines,
                [['C270106A', '1', 'L1_DAC'],
                 ['C060106F', '1', 'L1_HAC'],
                 ['C230998A-I3', '2', 'L23_BP'],
                 ['C020600C1', '4', 'L4_BTC'],
                 ['C060110A2', '5', 'L5_TTPC1'],
                 ['Fluo42_right', '6', 'L6_BPC']])


def test_duplicate_morphs():
    morph_db = morphdb.MorphologyDB()
    morph = MorphInfo(dict(name='morph0', mtype='L1_fake_type0', layer=1))
    ok_(morph_db.add_morph(morph))
    ok_(not morph_db.add_morph(morph))
    eq_(1, len(list(morph_db)))
