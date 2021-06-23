'''Utils'''
from functools import partial
from pathlib import Path
from typing import Optional

import pandas as pd
from deprecation import deprecated


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


@deprecated(details='Use `morph_tool.morphdb.MorphDB` instead.')
def neurondb_dataframe(neurondb: Path, morphology_dir: Optional[Path] = None) -> pd.DataFrame:
    '''Returns a DataFrame: [name, layer, mtype, use_axon, (optional) path]

    If read from an XML, additional columns maybe be present
    Args:
        filename: the neurondb.(dat|xml) file
        morphology_dir: (Optional) If passed, a column with the path to each morphology file
            will be added

    Raises: ValueError if the neurondb does not abide by the specification
    https://bbpteam.epfl.ch/documentation/projects/morphology-repair-workflow/latest/input_files.html#specification
    '''
    from morph_tool.morphdb import MorphDB  # pylint: disable=import-outside-toplevel,cyclic-import
    df = MorphDB.from_neurondb(neurondb, morphology_folder=morphology_dir).df
    columns = ['name', 'layer', 'mtype', 'use_axon']
    if not pd.isnull(df.path).all():
        columns.append('path')
    df = df[columns]
    return df
