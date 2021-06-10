"""
Caching morphology loader.
"""

import os

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

import morphio


def _ensure_startswith_point(file_ext):
    if file_ext.startswith("."):
        return file_ext
    else:
        return "." + file_ext


class MorphLoader(object):
    """
    Caching morphology loader.

    Args:
        base_dir: path to directory with morphology files
        file_ext: file extension to look for
        cache_size: size of LRU cache (if `None`, the cache can grow without bound,
            if 0, no lru_cache is used)
    """
    def __init__(self, base_dir, file_ext, cache_size=None):
        self.base_dir = base_dir
        self.file_ext = _ensure_startswith_point(file_ext)
        if cache_size == 0:
            self.get = self._get
        else:
            self.get = lru_cache(maxsize=cache_size)(self._get)

    def _get(self, name, options=None):
        filepath = os.path.join(self.base_dir, name + self.file_ext)
        if options is None:
            return morphio.Morphology(filepath)
        else:
            return morphio.Morphology(filepath, options)
