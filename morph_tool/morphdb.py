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
import collections
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd
import xmltodict
from more_itertools import always_iterable
from neurom.apps.morph_stats import extract_dataframe

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
                 **kwargs):
        '''MorphInfo ctor:

        Upon initialization, all repair attributes are set to True.

        Args:
            name: the morphology name (without the extension)
            mtype: the full mtype (ie 'L1_DAC:A')
            layer: the layer (accepts both string and int)
            kwargs: The following keyword arguments are also supported:

                - use_axon
                - use_dendrites
                - axon_repair
                - dendrite_repair
                - basal_dendrite_repair
                - tuft_dendrite_repair
                - oblique_dendrite_repair
                - unravel
                - use_for_stats
                - axon_inputs
                - path
                - label
        '''
        self.name = name
        self.mtype = mtype
        if self.MTYPE_SEPARATOR in mtype:
            self.mtype_no_subtype, self.msubtype = mtype.split(self.MTYPE_SEPARATOR)
        else:
            self.mtype_no_subtype, self.msubtype = self.mtype, ''
        self.layer = str(layer or '')

        for attr in BOOLEAN_REPAIR_ATTRS:
            setattr(self, attr, kwargs.get(attr, True))

        self.axon_inputs = kwargs.get('axon_inputs', [])
        self.path = kwargs.get('path')
        self.label = kwargs.get('label')

    @classmethod
    def _from_xmldict(cls, item: Dict[str, Any]):
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

        def is_true(repair: Dict[str, Any], key: str) -> bool:
            '''Parse a string representing a boolean repair flag and returns its boolean value
            Unless clearly stated as false, missing tags default to True
            - According to Eilif, an empty use_axon (corresponding to a null in the database)
              means that the axon is supposed to be used
            - dendrite_repair       defaults to True in BlueRepairSDK
            - basal_dendrite_repair defaults to True in BlueRepairSDK
            - unravel: well I guess we always want to do it
            '''
            el = repair.get(key)
            if el not in {'true', 'false', 'True', 'False', None}:
                raise ValueError(f'Invalid XML element {key} has invalid value: {el}\n'
                                 'Allowed values:\n'
                                 '- empty tag (which is equivalent to True)\n'
                                 '- true\n'
                                 '- True\n'
                                 '- false\n'
                                 '- False')
            return el in (None, 'true', 'True', )

        repair = item.get('repair', {})
        for attr in BOOLEAN_REPAIR_ATTRS:
            setattr(morph, attr, is_true(repair, attr))

        # "always_iterable" deals with <axoninput> not being interpreted
        # as a list if there is a single entry <axoninput> entry in the XML.
        morph.axon_inputs = list(always_iterable(
            (repair.get('axon_sources') or {}).get('axoninput', [])
        ))

        return morph

    @property
    def row(self) -> List:
        '''Flattened data structude ready to be used by a dataframe.'''
        return [getattr(self, attr) for attr in COLUMNS]

    def __repr__(self):
        return (f'MorphInfo(name={self.name!r}, mtype={self.mtype!r}, layer={self.layer!r}, '
                f'label={self.label!r})')


class MorphDB:
    '''A MorphInfo container.
    It takes care of maintaining uniqueness of the MorphInfo element
    and methods to write neurondb to various format (xml, dat, csv)
    '''

    def __init__(self, morph_info_seq: Optional[Iterable[MorphInfo]] = None):
        """Constructor of MorphDB.

        Args:
            morph_info_seq: an optional sequence of MorphInfo objects
        """
        self.df = pd.DataFrame([morph_info.row for morph_info in (morph_info_seq or ())],
                               columns=COLUMNS)
        MorphDB._sanitize_df_types(self.df)

    @classmethod
    def _from_neurondb_dat(cls, neurondb, morph_paths, label):
        '''Private constructor from neuronDB.dat files

        The equivalent public method is MorphDB.from_neurondb
        '''
        obj = cls()
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

    @classmethod
    def _from_neurondb_xml(cls, neurondb, morph_paths, label):
        obj = MorphDB()
        with neurondb.open() as fd:
            content = fd.read()
            neurondb = xmltodict.parse(content)

        morphologies = neurondb['neurondb']['listing']['morphology']

        # Case where there is a single <morphology></morphology> tag in the XML
        if not isinstance(morphologies, list):
            morphologies = [morphologies]
        morphologies = filter(None, morphologies)
        morphologies = list(map(MorphInfo._from_xmldict, morphologies))  # noqa, pylint: disable=protected-access

        for morph in morphologies:
            morph.label = label
            morph.path = morph_paths.get(morph.name)

        obj.df = MorphDB._create_dataframe(morphologies)
        MorphDB._sanitize_df_types(obj.df)
        return obj

    @classmethod
    def from_neurondb(cls,
                      neurondb: Union[Path, str],
                      label: str = 'default',
                      morphology_folder: Optional[Union[Path, str]] = None):
        '''Builds a MorphologyDB from a neurondb.(xml|dat) file
        Args:
            neurondb: path to a neurondb.(xml|dat) file
            label: a unique label to mark all morphologies coming from this neurondb
            morphology_folder: the location of the morphology files, if None it will default
                to the neurondb folder

        Raises: ValueError if the neurondb does not abide by the specification
        https://bbpteam.epfl.ch/documentation/projects/morphology-repair-workflow/latest/input_files.html#specification

        ..note:: missing keys are filled with `True` values
        '''
        neurondb = Path(neurondb)
        if morphology_folder:
            morphology_folder = Path(morphology_folder)
        else:
            morphology_folder = neurondb.parent.resolve()

        morph_paths = {path.stem: path for path in iter_morphology_files(morphology_folder)}

        if neurondb.suffix.lower() == '.dat':
            return cls._from_neurondb_dat(neurondb, morph_paths, label)
        else:
            return cls._from_neurondb_xml(neurondb, morph_paths, label)

    @classmethod
    def from_folder(cls,
                    morphology_folder: Union[Path, str],
                    mtypes: Iterable[Tuple[str, str]],
                    label: str = 'default',
                    extension: Optional[str] = None):
        '''Factory method to create a MorphDB object from a folder containing morphologies

        Args:
            morphology_folder: a folder containing morphologies
            mtypes: a sequence of 2-tuples (morphology name, mtype)
            label: (optional) a group label to be used to identify the morphlogies from this folder
            ext: Specify the morphology format to consider, if the folder contains multiple formats

        Raises: ValueError if the folder contains multiple files with the same name but
        different extensions and the extension argument has not been provided
        '''
        files = list(iter_morphology_files(Path(morphology_folder),
                                           extensions={extension} if extension else None))
        if not extension:
            duplicates = [item for item, count in
                          collections.Counter(path.stem for path in files).items()
                          if count > 1]
            if duplicates:
                raise ValueError(
                    f'Folder {morphology_folder} have multiple morphologies with the same '
                    'name but different extensions. This is not supported.\n'
                    f'Duplicate morphogies: {duplicates}\n\n'
                    'Please provide the extension to use with the arguement: extension')
        paths = {path.stem: path for path in files}
        return MorphDB(MorphInfo(name, mtype, label=label, path=paths[name])
                       for name, mtype in mtypes)

    def _to_xmldict(self):
        '''Transform the data to a xmldict compatible dictionary'''
        def get_repair_attr(morph, attr):
            if attr == 'axon_inputs':
                return 'axon_sources', {'axoninput': getattr(morph, attr)}
            else:
                return attr, getattr(morph, attr)

        return {'neurondb': {'listing': {
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

    def write(self, output_path: Union[Path, str]):
        '''Write the neurondb file to XML or DAT format

        Args:
            output_path: the output path
        '''
        output_path = Path(output_path)
        ext = output_path.suffix.lower()
        if ext == '.dat':
            self.df[['name', 'layer', 'mtype']].to_csv(output_path, sep=' ', header=False,
                                                       index=False)
        elif ext == '.xml':
            with output_path.open('w') as fd:
                fd.write(xmltodict.unparse(self._to_xmldict(), pretty=True))
        else:
            raise ValueError(f'Unsupported neurondb extensions ({ext}).'
                             ' Should be one of: (xml,csv,dat)')

    def features(self, config: Dict, n_workers=1):
        '''Returns a dataframe containing morphometrics and neurondb information

        Args:
            config: a NeuroM morph_stas config.
                See https://neurom.readthedocs.io/en/latest/morph_stats.html for more information
            n_workers: the number of workers to use to perform the computations

        Returns:
            A jointure dataframe between `neurom.stats.extract_dataframe` and `self.df`

        Raises:
            ValueError: if `self.df` has undefined (ie. None) paths
        '''
        missing_morphs = self.df[self.df.path.isnull()].name.to_list()
        if missing_morphs:
            raise ValueError(
                f'DataFrame has morphologies with undefined filepaths: {missing_morphs}')

        df = self.df.copy().reset_index(drop=True)
        df.columns = pd.MultiIndex.from_product((["properties"], df.columns.values))
        stats = extract_dataframe(df['properties', 'path'], config, n_workers)
        return df.join(stats.drop(columns='name', level=1), how='inner')

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
        obj += self
        obj += other
        return obj

    def __iadd__(self, other):
        if isinstance(other, MorphDB):
            self.df = pd.concat([self.df, other.df])
            MorphDB._sanitize_df_types(self.df)
        else:
            raise TypeError(f'Must be MorphDB or a sequence of MorphInfo, not {type(other)}')

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
        MorphDB._sanitize_df_types(df)
        return df

    @staticmethod
    def _sanitize_df_types(df):
        '''Set up the proper types for each columns

        Args:
            df: the dataframe to be sanitized
        '''
        df.astype({key: bool for key in BOOLEAN_REPAIR_ATTRS}, copy=False)
        df["axon_inputs"] = df["axon_inputs"].apply(tuple)
