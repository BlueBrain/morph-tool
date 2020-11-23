'''Utils'''
from functools import partial
from pathlib import Path
from typing import Optional

import pandas as pd
import xmltodict


def is_morphology(filename, extensions=None):
    '''Returns True if the extension is supported'''
    extensions = extensions or {'asc', 'h5', 'swc'}
    return Path(filename).suffix[1:].lower() in extensions


def iter_morphology_files(folder, recursive=False, extensions=None):
    '''Iterator that returns path to morphology files'''
    extensions = extensions or {'asc', 'h5', 'swc'}
    if recursive:
        files = Path(folder).rglob('*')
    else:
        files = Path(folder).glob('*')
    return filter(partial(is_morphology, extensions=extensions), files)


def find_morph(folder: Path, stem: str) -> Optional[Path]:
    '''Returns the path to a morphology in morphology_dir matching stem.

    If no morphology is found, returns None
    '''
    for ext in {'.asc', '.ASC', '.h5', '.H5', '.swc', '.SWC'}:
        path = folder / (stem + ext)
        if path.exists():
            return path
    return None


def _ensure_list(data):
    '''Returns the data if already a list else [data].

    xmltodict.parse will represent the data as a list if it sees
    the same tag multiple times in a row. For example,
    <morphology></morphology> will appear multiple times in neuronDB.xml

    However if the neurondb has a single <morphology></morphology> tag the data
    won't be represented as a list so this function enforce it

    Args:
        data: the python dict coming from xmltodict.parse.

    '''
    if isinstance(data, list):
        return data
    return [data]


def neurondb_dataframe(neurondb: Path, morphology_dir: Optional[Path] = None) -> pd.DataFrame:
    '''Returns a DataFrame: [name, layer, mtype]

    If read from an XML, additional columns maybe be present
    Args:
        filename: the neurondb.(dat|xml) file
        morphology_dir: (Optional) If passed, a column with the path to each morphology file
            will be added
    '''
    if neurondb.suffix.lower() == '.dat':
        columns = ['name', 'layer', 'mtype']
        df = pd.read_csv(neurondb, sep=r'\s+', names=columns, usecols=range(len(columns)))
        df.layer = df.layer.astype('str')
    elif neurondb.suffix.lower() == '.xml':
        with neurondb.open() as fd:
            neuronDB = xmltodict.parse(fd.read())

        rows = list()
        for morph in _ensure_list(neuronDB["neurondb"]["listing"]["morphology"]):
            use_axon = morph.get("repair", {}).get("use_axon")
            assert use_axon in {'true', 'false', 'True', 'False', None}
            rows.append(
                (morph["name"],
                 morph["layer"],
                 morph["mtype"] + (":" + morph["msubtype"] if morph["msubtype"] else ""),
                 # According to Eilif, an empty use_axon (corresponding to a null in the database)
                 # means that the axon is supposed to be used
                 (use_axon in {'true', 'True', None})))

        df = pd.DataFrame(data=rows, columns=['name', 'layer', 'mtype', 'use_axon'])
    else:
        raise ValueError(f'Unrecognized extension for neurondb file: {neurondb}')

    if morphology_dir:
        df['path'] = df.apply(lambda row: find_morph(morphology_dir, row['name']), axis=1)

    return df
