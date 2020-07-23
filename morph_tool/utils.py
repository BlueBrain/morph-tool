'''Utils'''
from functools import partial
from pathlib import Path

import pandas as pd
from morph_tool.morphdb import MorphologyDB


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


def neurondb_dataframe(filename: Path) -> pd.DataFrame:
    '''Returns a DataFrame: [name, layer, mtype]

    If read from an XML, additional columns maybe be present
    Args:
        filename: the neurondb.(dat|xml) file
    '''
    if filename.suffix.lower() == '.dat':
        columns = ['name', 'layer', 'mtype']
        df = pd.read_csv(filename, sep=r'\s+', names=columns, usecols=range(len(columns)))
        df.layer = df.layer.astype('str')
        return df

    if filename.suffix.lower() == '.xml':
        return MorphologyDB(filename).df

    raise ValueError(f'Unrecognized extension for neurondb file: {filename}')
