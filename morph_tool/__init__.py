'''Implementation of the morphology file converter'''
from morph_tool.version import __version__


class MorphToolException(Exception):
    '''MorphTool exception'''


class NoAxonException(MorphToolException):
    '''MorphTool exception'''


class NoDendriteException(MorphToolException):
    '''MorphTool exception'''


class TooManyAxonsException(MorphToolException):
    '''MorphTool exception'''
