'''This module provides the MorphIO diff functionality that can be used
to see if two morphologies are the same or not'''
import numpy as np
from morphio import Morphology


class DiffResult:
    '''
    An object that, when casted as a boolean, is equivalent to True when morphologies differ
    Additional information about why they differ is stored in DiffResult.info
    '''

    def __init__(self, is_different, info=None):
        '''The DiffResult ctor

        Args:
            is_different (bool): are the morphologies different
            info (str): Additional information about how morphologies differ
        '''
        self._is_different = is_different
        self.info = info

    def __bool__(self):
        '''Returns True if morphologies differ'''
        return self._is_different

    def __nonzero__(self):
        '''Returns True if morphologies differ, but for python 2'''
        return self.__bool__()


def diff(morph1, morph2, rtol=1.e-5, atol=1.e-8):
    '''
    Returns a DiffResult object that is equivalent to True when morphologies differ
    Additional information about why they differ is stored in DiffResult.info

    Morphologies with different formats can be compared.

    Morphologies are considered different if one of the following property differ:
    - number of root sections
    - sections type
    - sections point array
    - sections diameter array
    - sections perimeter array
    - sections number of children

    The soma are NOT taken into consideration

    Args:
        morph1 (str|morphio.Morphology|morphio.mut.Morphology): a morphology
        morph2 (str|morphio.Morphology|morphio.mut.Morphology): a morphology
        rtol (float): the relative tolerance used for comparing points (see numpy.isclose help)
        atol (float): the absolute tolerance used for comparing points (see numpy.isclose help)
    '''

    if not isinstance(morph1, Morphology):
        morph1 = Morphology(morph1)
    if not isinstance(morph2, Morphology):
        morph2 = Morphology(morph2)

    if len(morph1.root_sections) != len(morph2.root_sections):
        return DiffResult(True,
                          'Both morphologies have a different number of root sections')

    for section1, section2 in zip(morph1.iter(), morph2.iter()):
        for attrib in ['points', 'diameters', 'perimeters']:
            val1, val2 = getattr(section1, attrib), getattr(section2, attrib)
            if val1.shape != val2.shape:
                return DiffResult(True,
                                  'Attributes Section.{} of:\n'
                                  '{}\n{}\nhave different shapes: {} vs {}'.format(
                                      attrib, section1, section2, val1.shape, val2.shape))
            is_close = np.isclose(val1, val2, rtol=rtol, atol=atol)
            if not is_close.all():
                first_diff_index = np.where(~is_close)[0][0]

                return DiffResult(True,
                                  f'Attributes Section.{attrib} of:\n'
                                  f'{section1}\n{section2}\nhave the same shape '
                                  'but different values\n'
                                  f'Vector {attrib} differs at index {first_diff_index}: '
                                  f'{val1[first_diff_index]} != {val2[first_diff_index]}')
        if section1.type != section2.type:
            return DiffResult(True, '{} and {} have different section types'.format(
                section1, section2))
        if len(section1.children) != len(section2.children):
            return DiffResult(True,
                              '{} and {} have a different number of children'.format(
                                  section1, section2))

    return DiffResult(False)
