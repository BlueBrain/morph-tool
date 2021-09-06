"""Constants module.

Constants for names that are required to be in columns of ``synapses`` arg of
``morph_tool.plot`` package. These constants are the same as edge properties of `Sonata
<https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md#
edges---optional-reserved-attributes>`__. So you can use these constants, their values as plain
strings or constants from `bluepysnap <https://github.com/BlueBrain/snap/>`__ library:

An example of using constants or plain strings:
    .. literalinclude:: ../../../examples/dendrogram.py
        :pyobject: plain_example

An example of using with a circuit and bluepysnap constants:
    .. literalinclude:: ../../../examples/dendrogram.py
        :pyobject: circuit_example

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
