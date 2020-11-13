'''Implementation of the morphology file converter'''
from morphio import LogLevel, Morphology

from morph_tool.apical_point import apical_point_position, apical_point_section_segment
from morph_tool.converter import convert
from morph_tool.morphio_diff import diff
from morph_tool.spatial import point_to_section_segment


class MorphToolException(Exception):
    '''MorphTool exception'''


class NoAxonException(MorphToolException):
    '''MorphTool exception'''


class NoDendriteException(MorphToolException):
    '''MorphTool exception'''
