'''Module providing utility functions for the tests'''
import shutil
import tempfile
from contextlib import contextmanager


@contextmanager
def setup_tempdir(prefix, no_cleanup=False):
    '''Context manager returning a temporary directory'''
    temp_dir = tempfile.mkdtemp(prefix=prefix)
    try:
        yield temp_dir
    finally:
        if not no_cleanup:
            shutil.rmtree(temp_dir)
