'''Implementation of the morphology file converter'''
from morphio import LogLevel, Morphology
from morphio import diff as diff_

from morph_tool.version import __version__


def diff(morph1, morph2, verbose_level=LogLevel.info):
    '''
    Returns True if MORPH1 and MORPH2 differ, False otherwise.

    Morphologies with different formats can be compared.

    Args:
        morph1 (str|morphio.Morphology): a morphology
        morph2 (str|morphio.Morphology): a morphology
        verbose_level (morphio.LogLevel): a verbosity level
    '''
    if not isinstance(morph1, Morphology):
        morph1 = Morphology(morph1)

    if not isinstance(morph2, Morphology):
        morph2 = Morphology(morph2)

    return diff_(morph1, morph2, verbose_level)


class MorphToolException(Exception):
    '''MorphTool exception'''


class NoAxonException(MorphToolException):
    '''MorphTool exception'''


class NoDendriteException(MorphToolException):
    '''MorphTool exception'''


class TooManyAxonsException(MorphToolException):
    '''MorphTool exception'''
