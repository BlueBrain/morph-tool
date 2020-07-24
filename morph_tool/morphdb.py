'''Provides two classes to work with neurondb files

The main use is to provide a neuronDB loader and a method to retrieve
information as a dataframe:

MorphologyDB('neurondb.xml').df
'''
import json
import logging
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional

import pandas as pd
import xmltodict

L = logging.getLogger(__name__)

NEURONDB_XML = 'neuronDB.xml'
MTYPE_MSUBTYPE_SEPARATOR = ':'


class MorphInfo:
    '''A class the contains information about a morphology.

    Its role is to abstract away the raw data.
    '''

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

    COLUMNS = ['name', 'mtype', 'msubtype', 'fullmtype',
               'layer'] + BOOLEAN_REPAIR_ATTRS + ['axon_inputs']

    def __init__(self, item: Dict[str, Any]):
        '''MorphInfo ctor.

        Args:
            item: A dictionnary that represents the content of the XML file.
                  The only mandatory keys are [name, mtype, layer]
        '''
        self.item = item
        self.name = item['name']
        self.mtype = item['mtype']
        self.layer = str(item['layer'])
        self.msubtype = item.get('msubtype') or ''
        self.fullmtype = MTYPE_MSUBTYPE_SEPARATOR.join(filter(None, [self.mtype, self.msubtype]))

        def is_true(el: Optional[str]) -> bool:
            '''Parse a string representing a boolean repair flag and returns its boolean value

            Unless clearly stated as false, missing tags default to True

            - According to Eilif, an empty use_axon (corresponding to a null in the database)
              means that the axon is supposed to be used

            - dendrite_repair       defaults to True in BlueRepairSDK
            - basal_dendrite_repair defaults to True in BlueRepairSDK
            - unravel: well I guess we always want to do it
            '''
            return el in (None, '', 'true', 'True', )

        repair = item.get('repair', {})
        for attr in MorphInfo.BOOLEAN_REPAIR_ATTRS:
            setattr(self, attr, is_true(repair.get(attr)))

        # Case where there is a single <axoninput></axoninput> tag in the XML
        self.axon_inputs = repair.get('axon_sources', {}).get('axoninput', [])
        if not isinstance(self.axon_inputs, list):
            self.axon_inputs = [self.axon_inputs]

        # lineage information
        self.dendrite_donor = item.get('parent') or item.get('dendrite')
        self.axon_donor = item.get('parent') or item.get('axon')

    def __hash__(self):
        return hash((self.name, self.mtype, self.layer))

    @property
    def data(self) -> Dict:
        '''Data that matter to generate the neurondb.xml'''
        return {
            'name': self.name,
            'mtype': self.mtype,
            'msubtype': self.msubtype,
            'layer': self.layer,
            'repair': {attr: getattr(self, attr)
                       for attr in MorphInfo.BOOLEAN_REPAIR_ATTRS + ['axon_inputs']}
        }

    @property
    def row(self) -> List:
        '''Flattened data structude ready to be used by a dataframe.'''
        return [getattr(self, attr) for attr in MorphInfo.COLUMNS]

    def __repr__(self):
        return f'MorphInfo(name={self.name}, mtype={self.mtype}, layer={self.layer})'


class MorphologyDB(object):
    '''A MorphInfo container.

    It takes care of maintaining unicity of the MorphInfo element
    and methods to write neurondb to various format (xml, dat, csv)
    '''

    def __init__(self,
                 neurondb: Path = None,
                 morph_info_filter: Callable[[MorphInfo], bool] = None):
        '''Builds a MorphologyDB from a neurondb.xml file

        Args:
            neurondb: path to a neurondb.xml file
            morph_info_filter: a filter function to be applied to each MorphInfo element
        '''
        self.lineage = {}

        if neurondb is not None:
            with open(neurondb) as fd:
                neurondb = xmltodict.parse(fd.read())

            morphologies = neurondb['neurondb']['listing']['morphology']

            # Case where there is a single <morphology></morphology> tag in the XML
            if not isinstance(morphologies, list):
                morphologies = [morphologies]
            morphologies = filter(None, morphologies)
            morphologies = map(MorphInfo, morphologies)
            if morph_info_filter is not None:
                morphologies = filter(morph_info_filter, morphologies)
        else:
            morphologies = []

        self.morphologies = list(morphologies)
        self._known_morphologies = set(map(hash, self.morphologies))

    def __iter__(self):
        return iter(self.morphologies)

    def remove_morphs(self, removed_morphs: Iterable[MorphInfo]) -> None:
        '''Removes morphologies from the database.'''
        removed_morphs = set(removed_morphs)
        self.morphologies = [morph for morph in self.morphologies
                             if morph.name not in removed_morphs]

        self._known_morphologies = set(map(hash, self.morphologies))

        for name in removed_morphs:
            self.lineage.pop(name, None)

    def add_morph(self, morph_info: MorphInfo) -> bool:
        '''Add a morphology to the database.

        If the name, mtype, layer combination already exists, then nothing is added

        Return:
            True if morphology added, False otherwise
        '''
        key = hash(morph_info)
        if key in self._known_morphologies:
            return False
        self._known_morphologies.add(key)
        self.morphologies.append(morph_info)

        lineage = {'dendrite': morph_info.dendrite_donor, 'axon': morph_info.axon_donor}
        lineage = {k: v for k, v in lineage.items() if v is not None}
        if lineage:
            self.lineage[morph_info.name] = lineage
        return True

    @property
    def df(self):
        '''Returns a pandas.DataFrame view of the data with the following columns:

           'name': the morpho name (without extension)
           'mtype': the mtype (without msubtype),
           'msubtype': the msubtype (the part of the fullmtype after the ":")
           'fullmtype': the full mtype (mtype:msubtype)
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
           'use_for_stats': ???
           'axon_inputs': the list of morphologies whose axon can be grafted on this morphology
        '''
        return pd.DataFrame([morph.row for morph in self.morphologies],
                            columns=MorphInfo.COLUMNS)

    @property
    def data(self):
        '''The neurondb relevant data'''
        return {'neurondb':
                {'overview': {'count': len(self.morphologies)},
                 'listing': {'morphology': [morph.data for morph in self.morphologies]}}}

    def write(self,
              output_path: Path,
              filt: Dict[str, bool] = None):
        '''Write the neurondb file to XML, DAT or CSV format'''
        output_path = Path(output_path)
        ext = output_path.suffix.lower()[1:]
        if ext == 'csv':
            self._write_neurondb(output_path, sep=',', write_header=True, filt=filt)
        elif ext == 'dat':
            self._write_neurondb(output_path, sep=' ', write_header=False, filt=filt)
        elif ext == 'xml':
            with output_path.open('w') as fd:
                fd.write(xmltodict.unparse(self.data, pretty=True))
        else:
            raise ValueError(f'Unsupported neurondb extensions ({ext}).'
                             ' Should be one of: (xml,csv,dat)')
        return output_path

    def _write_neurondb(self, output_path: Path, sep=' ', write_header=False,
                        filt: Dict[str, bool] = None):
        if filt is None:
            filt = {}
        with open(output_path, 'w') as fd:
            if write_header:
                fd.write(sep.join(['name', 'layer', 'mtype']) + "\n")
            for morph in self:
                if ('use_axon' in filt) and (morph.use_axon != filt['use_axon']):
                    continue
                if ('use_dendrites' in filt) and (morph.use_dendrites != filt['use_dendrites']):
                    continue
                fd.write(sep.join([morph.name, morph.layer, morph.mtype]) + "\n")

    def write_lineage(self, output_path: Path, name='lineage.json'):
        '''write the lineage to the output_path/name'''
        full_path = Path(output_path, name)
        with open(full_path, 'w') as fd:
            json.dump(self.lineage, fd, indent=2)
        return full_path
