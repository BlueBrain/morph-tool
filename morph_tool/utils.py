'''Utils'''
from functools import partial
from pathlib import Path

import pandas as pd


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

    Args:
        filename: the neurondb.dat file
    '''
    df = pd.read_csv(filename, sep=' ', names=['name', 'layer', 'mtype'])
    df.layer = df.layer.astype('str')
    return df
