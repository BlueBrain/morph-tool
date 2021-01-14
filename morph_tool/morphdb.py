'''Provides two classes to work with neurondb files
The main use is to provide a neuronDB loader and a method to retrieve
information as a dataframe.

Example:

old_release = MorphDB.from_neurondb(Path("/realease/one/neuronDB.xml"), label='old-release')
new_release = MorphDB.from_neurondb(Path("/realease/two/neuronDB.xml"), label='new-release')
total = old_release + new_release

# neurondb information
print(total.df)

# combined information (neurondb + morph_stats)
print(total.features({'neurite': {'section_lengths': ['max']}}))
'''
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import xmltodict
from neurom.apps.morph_stats import extract_dataframe
from toolz import isiterable

from morph_tool.utils import iter_morphology_files

L = logging.getLogger(__name__)


# Attributes found in the <repair></repair> block of the XML
BOOLEAN_REPAIR_ATTRS = [
    'use_axon',
    'use_dendrites',
    'axon_repair',
    'dendrite_repair',
    'basal_dendrite_repair',
    'tuft_dendrite_repair',
    'oblique_dendrite_repair',
    'unravel',
    'use_for_stats',
]

COLUMNS = ['name', 'mtype', 'msubtype', 'mtype_no_subtype',
           'layer', 'label', 'path'] + BOOLEAN_REPAIR_ATTRS + ['axon_inputs']


class MorphInfo:
    '''A class the contains information about a morphology.
    Its role is to abstract away the raw data.
    '''

    MTYPE_SEPARATOR = ':'

    def __init__(self, name: str, mtype: str, layer: Optional[Union[str, int]] = None,
                 label: str = None):
        '''MorphInfo ctor:

        Upon initialization, all repair attributes are set to True.

        Args:
            name: the morphology name (without the extension)
            mtype: the full mtype (ie 'L1_DAC:A')
            layer: the layer (accepts both string and int)
            label: (optional) a group label to be used to identify multiple groups of morphologies
                added to the same collection.
        '''
        self.name = name
        self.mtype = mtype
        if self.MTYPE_SEPARATOR in mtype:
            self.mtype_no_subtype, self.msubtype = mtype.split(self.MTYPE_SEPARATOR)
        else:
            self.mtype_no_subtype, self.msubtype = self.mtype, ''
        self.layer = str(layer or '')

        for attr in BOOLEAN_REPAIR_ATTRS:
            setattr(self, attr, True)

        self.axon_inputs = []
        self.path = None
        self.label = label

    @classmethod
    def _from_xml(cls, item: Dict[str, Any]):
        '''MorphInfo ctor.
        Args:
            item: A dictionary that represents the content of the XML file.
                  The only mandatory keys are [name, mtype, layer]
        '''
        morph = cls(
            item['name'],
            cls.MTYPE_SEPARATOR.join(filter(None, [item['mtype'], item.get('msubtype')])),
            item['layer']
        )

        def is_true(el: Optional[str]) -> bool:
            '''Parse a string representing a boolean repair flag and returns its boolean value
            Unless clearly stated as false, missing tags default to True
            - According to Eilif, an empty use_axon (corresponding to a null in the database)
              means that the axon is supposed to be used
            - dendrite_repair       defaults to True in BlueRepairSDK
            - basal_dendrite_repair defaults to True in BlueRepairSDK
            - unravel: well I guess we always want to do it
            '''
            assert el in {'true', 'false', 'True', 'False', None}, f'Invalid element: {el}'
            return el in (None, '', 'true', 'True', )

        repair = item.get('repair', {})
        for attr in BOOLEAN_REPAIR_ATTRS:
            setattr(morph, attr, is_true(repair.get(attr)))

        axon_sources = repair.get('axon_sources')
        morph.axon_inputs = axon_sources.get('axoninput', []) if axon_sources else []

        # Case where there is a single <axoninput></axoninput> tag in the XML
        if not isinstance(morph.axon_inputs, list):
            morph.axon_inputs = [morph.axon_inputs]

        return morph

    @property
    def row(self) -> List:
        '''Flattened data structude ready to be used by a dataframe.'''
        return [getattr(self, attr) for attr in COLUMNS]

    def __repr__(self):
        return f'MorphInfo(name={self.name}, mtype={self.mtype}, layer={self.layer})'


class MorphDB(object):
    '''A MorphInfo container.
    It takes care of maintaining unicity of the MorphInfo element
    and methods to write neurondb to various format (xml, dat, csv)
    '''

    def __init__(self):
        self.df = pd.DataFrame(columns=COLUMNS)
        self.df.astype({key: bool for key in BOOLEAN_REPAIR_ATTRS}, copy=False)

    @classmethod
    def from_neurondb(cls,
                      neurondb: Path = None,
                      label: str = 'default',
                      morphology_folder: Optional[Path] = None):
        '''Builds a MorphologyDB from a neurondb.(xml|dat) file
        Args:
            neurondb: path to a neurondb.(xml|dat) file
            label: a unique label to mark all morphologies coming from this neurondb
            morphology_folder: the location of the morphology files, if None it will default
                to the neurondb folder

        ..note:: missing keys are filled with `True` values
        '''
        obj = MorphDB()

        if not morphology_folder:
            morphology_folder = neurondb.parent.resolve()

        morph_paths = {path.stem: path for path in iter_morphology_files(morphology_folder)}

        if neurondb.suffix.lower() == '.dat':
            columns = ['name', 'layer', 'mtype']
            obj.df = pd.read_csv(neurondb, sep=r'\s+', names=columns, usecols=range(len(columns)))

            fulltypes = obj.df.mtype.str.split(':', n=1, expand=True).fillna('')
            if len(fulltypes.columns) > 1:
                obj.df['mtype_no_subtype'] = fulltypes[0]
                obj.df['msubtype'] = fulltypes[1]
            else:
                obj.df['mtype_no_subtype'] = obj.df.mtype
                obj.df['msubtype'] = ''
            obj.df['label'] = label
            for missing_col in set(COLUMNS) - set(obj.df.columns):
                obj.df[missing_col] = None
            obj.df.layer = obj.df.layer.astype('str')
            obj.df['path'] = obj.df.name.map(morph_paths)
            obj.df = obj.df.reindex(columns=COLUMNS)
            for key in BOOLEAN_REPAIR_ATTRS:
                obj.df[key] = True
            obj.df['axon_inputs'] = [[] for _ in range(len(obj.df))]

            return obj

        with neurondb.open() as fd:
            content = fd.read()
            neurondb = xmltodict.parse(content)

        morphologies = neurondb['neurondb']['listing']['morphology']

        # Case where there is a single <morphology></morphology> tag in the XML
        if not isinstance(morphologies, list):
            morphologies = [morphologies]
        morphologies = filter(None, morphologies)
        morphologies = list(map(MorphInfo._from_xml, morphologies))  # noqa, pylint: disable=protected-access

        for morph in morphologies:
            morph.label = label
            morph.path = morph_paths.get(morph.name)

        obj.df = MorphDB._create_dataframe(morphologies)
        obj.df = obj.df.astype({key: bool for key in BOOLEAN_REPAIR_ATTRS})
        return obj

    @classmethod
    def from_folder(cls,
                    morphology_folder: Path,
                    mtypes: Dict[str, str],
                    label: str = 'default'):
        '''Factory method to create a MorphDB object from a folder containing morphologies

        Args:
            morphology_folder: a folder containing morphologies
            mtype: a dictionary of where key is a morphology name and value its mtype
            label: (optional) a group label to be used to identify the morphlogies from this folder
        '''

        obj = cls()
        obj += (MorphInfo(path.stem, mtypes[path.stem], label=label)
                for path in iter_morphology_files(morphology_folder))
        return obj

    def write(self, output_path: Path):
        '''Write the neurondb file to XML or DAT format'''
        ext = output_path.suffix.lower()
        if ext == '.dat':
            self.df[['name', 'layer', 'mtype']].to_csv(output_path, sep=' ', header=False,
                                                       index=False)
        elif ext == '.xml':
            def get_repair_attr(morph, attr):
                if attr == 'axon_inputs':
                    return 'axon_sources', {'axoninput': getattr(morph, attr)}
                else:
                    return attr, getattr(morph, attr)
            data = {'neurondb':
                    {'listing': {
                        'morphology': [{
                            'name': morph.name,
                            'mtype': morph.mtype_no_subtype,
                            'msubtype': morph.msubtype,
                            'layer': morph.layer,
                            'repair': dict(
                                get_repair_attr(morph, attr)
                                for attr in BOOLEAN_REPAIR_ATTRS + ['axon_inputs']
                            )
                        } for morph in self.df.itertuples()]}
                     }}
            with output_path.open('w') as fd:
                fd.write(xmltodict.unparse(data, pretty=True))
        else:
            raise ValueError(f'Unsupported neurondb extensions ({ext}).'
                             ' Should be one of: (xml,csv,dat)')
        return output_path

    def features(self, config: Dict, n_workers=1):
        '''Returns a dataframe containing morphometrics and neurondb information

        Args:
            config: a NeuroM morph_stas config.
                See https://neurom.readthedocs.io/en/latest/morph_stats.html for more information

        Returns:
            A jointure dataframe between `neurom.stats.extract_dataframe` and `self.df`

        Raises:
            ValueError: if `self.df` has undefined (ie. None) paths
        '''
        missing_morphs = self.df[self.df.path.isnull()].name.values
        if missing_morphs:
            raise ValueError(
                f'DataFrame has morphologies with undefined filepaths: {missing_morphs}')

        df = self.df.copy().reset_index()
        df.columns = pd.MultiIndex.from_product((["neuron"], df.columns.values))
        paths = df['neuron', 'path']

        stats = extract_dataframe(paths, config, n_workers).drop(columns='name', level=1)
        return df.join(stats, how='inner').drop(columns=('neuron', 'index'))

    def check_files_exist(self):
        '''Raises if `self.df.path` has None values or non existing paths'''
        missing_morphs = self.df[self.df.path.isnull()].name.values
        if missing_morphs:
            raise ValueError(
                f'DataFrame has morphologies with undefined filepaths: {missing_morphs}')

        for path in self.df.path:
            if not path.exists():
                raise ValueError(f'Non existing path: {path}')

    def __add__(self, other):
        obj = MorphDB()
        if isinstance(other, MorphDB):
            obj.df = pd.concat([self.df, other.df])
        elif isiterable(other):
            seq = list(other)
            if seq and isinstance(seq[0], MorphInfo):
                df = pd.DataFrame([morph_info.row for morph_info in seq], columns=COLUMNS)
                obj.df = pd.concat([self.df, df])
        else:
            raise TypeError(f'Must be MorphDB or a sequence of MorphInfo, not {type(other)}')

        obj.df = obj.df.astype({key: bool for key in BOOLEAN_REPAIR_ATTRS})
        return obj

    def __iadd__(self, other):
        if isinstance(other, MorphDB):
            self.df = pd.concat([self.df, other.df])
        elif isiterable(other):
            seq = list(other)
            if seq and isinstance(seq[0], MorphInfo):
                df = pd.DataFrame([morph_info.row for morph_info in seq], columns=COLUMNS)
                self.df = pd.concat([self.df, df])
        else:
            raise TypeError(f'Must be MorphDB or a sequence of MorphInfo, not {type(other)}')

        self.df = self.df.astype({key: bool for key in BOOLEAN_REPAIR_ATTRS})
        return self

    @staticmethod
    def _create_dataframe(morphologies: List[MorphInfo]):
        '''Returns a pandas.DataFrame view of the data with the following columns:
           'name': the morpho name (without extension)
           'mtype': the mtype with its subtype
           'msubtype': the msubtype (the part of the mtype after the ":")
           'mtype_no_subtype': the mtype without the msubtype
           'layer': the layer (as a string)
           # repair related columns
           'use_axon': states that the morphology's axon can be used as a donor for axon grafting
           'use_dendrites': states that this morphology can be used as a recipient for axon grafting
           'axon_repair': flag to activate axon repair
           'dendrite_repair': flag to activate dendrites repair
           'basal_dendrite_repair': flag to activate basal dendrites repair (dendrite_repair
                                    must be true as well)
           'tuft_dendrite_repair': flag to activate tuft dendrites repair (dendrite_repair
                                    must be true as well)
           'oblique_dendrite_repair': flag to activate oblique dendrites repair (dendrite_repair
                                    must be true as well)
           'unravel': flag to activate unravelling
           'axon_inputs': the list of morphologies whose axon can be grafted on this morphology
           # deprecated
           'use_for_stats': Legacy flag that was used to determine if an axon was suitable to be
                            used as a axoninput
        '''
        df = pd.DataFrame([morph.row for morph in morphologies], columns=COLUMNS)
        df.astype({key: bool for key in BOOLEAN_REPAIR_ATTRS}, copy=False)
        return df
