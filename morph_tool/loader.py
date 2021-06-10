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


def get_incomplete_sections(morph_path, tag="Incomplete"):
    """Read certain tags in ascii files (Incomplete/Normal) corresponding to leaves.

    WARNING: this function is fairly experimental, and not tested on all possible cases

    Args:
        morph_path (str): path to ascii morphology file
        tag (str): tag to detect

    Returns:
        list: list of morphio sections ids with this tag
    """
    from morph_tool.spatial import point_to_section_segment

    if not str(morph_path).lower().endswith('.asc'):
        raise Exception("This only works on .asc files")

    neuron = morphio.Morphology(morph_path)

    incomplete_sections = []
    with open(morph_path, "rb") as morph:
        _last_point = None

        for i, line in enumerate(morph):
            # skip badly encoded lines (not part of the morphology)
            try:
                _line = bytearray(line).decode()
            except UnicodeDecodeError:
                continue

            # update last red point
            if "(" in _line and ")" in _line:
                _last_point = _line

            # if tag is found on the line, record the section id of that point
            if tag in _line:
                incomplete_sections.append(
                    point_to_section_segment(
                        neuron,
                        [
                            float(_point)
                            for _point in _last_point.split("(")[1].split(")")[0].split(" ")
                            if _point != ""
                        ],
                    )[0]
                )

    return incomplete_sections
