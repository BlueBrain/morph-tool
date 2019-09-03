'''Implementation of the morphology file converter'''
from morphio import LogLevel, Morphology

from morph_tool.apical_point import get_apical_point
from morph_tool.morphio_diff import diff
from morph_tool.version import __version__


class MorphToolException(Exception):
    '''MorphTool exception'''


class NoAxonException(MorphToolException):
    '''MorphTool exception'''


class NoDendriteException(MorphToolException):
    '''MorphTool exception'''


class TooManyAxonsException(MorphToolException):
    '''MorphTool exception'''
