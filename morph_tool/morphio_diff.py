"""Functionality that can be used to see if two morphologies are the same or not."""
import numpy as np
from morphio import Morphology


class DiffResult:
    """An object that, when casted as a boolean, is equivalent to True when morphologies differ.

    Additional information about why they differ is stored in DiffResult.info.
    """

    def __init__(self, is_different, info=None):
        """The DiffResult constructor.

        Args:
            is_different (bool): are the morphologies different
            info (str): Additional information about how morphologies differ
        """
        self._is_different = is_different
        self.info = info

    def __bool__(self):
        """Returns True if morphologies differ."""
        return self._is_different

    def __nonzero__(self):
        """Returns True if morphologies differ, but for python 2."""
        return self.__bool__()


def diff(morph1, morph2, rtol=1.e-5, atol=1.e-8, *, skip_perimeters=False, all_diffs=False):
    """Returns a DiffResult object that is equivalent to True when morphologies differ.

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
        skip_perimeters (bool): do not check the perimeters if set to True
        all_diffs (bool): return all differences if set to True
    """
    # pylint: disable=too-many-branches
    # pylint: disable=too-many-locals
    # pylint: disable=too-many-return-statements
    if not isinstance(morph1, Morphology):
        morph1 = Morphology(morph1)
    if not isinstance(morph2, Morphology):
        morph2 = Morphology(morph2)

    diffs = []
    if len(morph1.root_sections) != len(morph2.root_sections):
        return DiffResult(True,
                          'Both morphologies have a different number of root sections')

    attr_list = ['points', 'diameters']
    if not skip_perimeters:
        attr_list.append('perimeters')

    for section1, section2 in zip(morph1.iter(), morph2.iter()):
        for attrib in attr_list:
            val1, val2 = getattr(section1, attrib), getattr(section2, attrib)
            if val1.shape != val2.shape:
                current_diff = (
                    f'Attributes Section.{attrib} of:\n'
                    f'{section1}\n'
                    f'{section2}\n'
                    f'have different shapes: {val1.shape} vs {val2.shape}'
                )
                if all_diffs:
                    diffs.append(current_diff)
                else:
                    return DiffResult(True, current_diff)
            is_close = np.isclose(val1, val2, rtol=rtol, atol=atol)
            if not is_close.all():
                first_diff_index = np.where(~is_close)[0][0]
                current_diff = (
                    f'Attributes Section.{attrib} of:\n'
                    f'{section1}\n{section2}\nhave the same shape '
                    'but different values\n'
                    f'Vector {attrib} differs at index {first_diff_index}: '
                    f'{val1[first_diff_index]} != {val2[first_diff_index]}'
                )
                if all_diffs:
                    diffs.append(current_diff)
                else:
                    return DiffResult(True, current_diff)
        if section1.type != section2.type:
            current_diff = f'{section1} and {section2} have different section types'
            if all_diffs:
                diffs.append(current_diff)
            else:
                return DiffResult(True, current_diff)
        if len(section1.children) != len(section2.children):
            current_diff = f'{section1} and {section2} have a different number of children'
            if all_diffs:
                diffs.append(current_diff)
            else:
                return DiffResult(True, current_diff)

    if all_diffs and len(diffs) > 0:
        return DiffResult(True, "\n\n".join(diffs))
    return DiffResult(False)
