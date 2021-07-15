"""Constants for names that are required to be in columns of ``synapses`` arg of
``morph_tool.plot`` package.

You either use these constants for columns names of ``synapses`` arg or you can use their
values directly as plain strings.

An example from ``examples`` directory::
    .. literalinclude:: ../../../examples/dendrogram.py
        :pyobject: plain_example

"""

#: Contains string `@target_node` that must be used as column name in `synapses` arg.
TARGET_NODE_ID = '@target_node'
#: Contains string `@source_node` that must be used as column name in `synapses` arg.
SOURCE_NODE_ID = '@source_node'

#: Contains string `afferent_section_id` that must be used as column name in `synapses` arg.
POST_SECTION_ID = 'afferent_section_id'
#: Contains string `afferent_section_pos` that must be used as column name in `synapses` arg.
POST_SECTION_POS = 'afferent_section_pos'
#: Contains string `afferent_segment_id` that must be used as column name in `synapses` arg.
POST_SEGMENT_ID = 'afferent_segment_id'
#: Contains string `afferent_segment_offset` that must be used as column name in `synapses` arg.
POST_SEGMENT_OFFSET = 'afferent_segment_offset'

#: Contains string `efferent_section_id` that must be used as column name in `synapses` arg.
PRE_SECTION_ID = 'efferent_section_id'
#: Contains string `efferent_section_pos` that must be used as column name in `synapses` arg.
PRE_SECTION_POS = 'efferent_section_pos'
#: Contains string `efferent_segment_id` that must be used as column name in `synapses` arg.
PRE_SEGMENT_ID = 'efferent_segment_id'
#: Contains string `efferent_segment_offset` that must be used as column name in `synapses` arg.
PRE_SEGMENT_OFFSET = 'efferent_segment_offset'
