'''Utils'''
import os
from functools import partial


def is_morphology(filename, extensions=None):
    '''Returns True if the extension is supported'''
    s = filename.split('.')
    if len(s) < 2:
        return False
    extensions = extensions or {'asc', 'h5', 'swc'}
    return s[-1].lower() in extensions


def iter_morphology_files(folder, recursive=False, extensions=None):
    '''Iterator that returns path to morphology files'''
    extensions = extensions or {'asc', 'h5', 'swc'}
    if recursive:
        files = (os.path.join(root, file) for root, _, files in os.walk(folder)
                 for file in files)
    else:
        files = (os.path.join(folder, f) for f in os.listdir(folder))
    return filter(partial(is_morphology, extensions=extensions), files)
