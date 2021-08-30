import numpy as np
import pandas as pd
import neurom as nm
from morph_tool.plot import morphology
from morph_tool.plot import consts


def example_afferent():
    m = nm.load_morphology('dendrogram_plain_example.swc')
    data = np.array([
        [4, 0, 0.81408846, 'additional value'],
        [6, 2, 0.545983203, 'additional value'],
        [1, 0, 0.4290702, 'additional value'],
    ])
    columns = [consts.POST_SECTION_ID, consts.POST_SEGMENT_ID, consts.POST_SEGMENT_OFFSET,
               'additional_data']
    synapses = pd.DataFrame(data, columns=columns)

    fig = morphology.draw(m, synapses)
    fig.show()


def example_efferent():
    m = nm.load_morphology('dendrogram_plain_example.swc')
    data = np.array([
        [4, 0, 0.81408846, 'additional value'],
        [6, 2, 0.645983203, 'additional value'],
        [1, 0, 0.4290702, 'additional value'],
    ])
    columns = [consts.PRE_SECTION_ID, consts.PRE_SEGMENT_ID, consts.PRE_SEGMENT_OFFSET,
               'additional_data']
    synapses = pd.DataFrame(data, columns=columns)

    fig = morphology.draw(m, synapses)
    fig.show()


def example_afferent_efferent():
    m = nm.load_morphology('dendrogram_plain_example.swc')
    data = np.array([
        [4, 0, 0.81408846, np.nan, np.nan, np.nan],
        [6, 2, 0.145983203, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, 1, 0, 0.4290702],
        [np.nan, np.nan, np.nan, 1, 0, 0.29180855],
        [np.nan, np.nan, np.nan, 1, 0, 0.68261607],
    ])
    columns = [consts.POST_SECTION_ID, consts.POST_SEGMENT_ID, consts.POST_SEGMENT_OFFSET,
               consts.PRE_SECTION_ID, consts.PRE_SEGMENT_ID, consts.PRE_SEGMENT_OFFSET]
    synapses = pd.DataFrame(data, columns=columns)

    fig = morphology.draw(m, synapses)
    fig.show()


def circuit_example():
    """Example that shows how to plot a morphology with synapses from a Sonata circuit.

    To make this example work, you would need a proper SONATA circuit.
    """
    from bluepysnap import Circuit
    from bluepysnap.bbp import Synapse
    circuit = Circuit('/path/to/sonata_circuit_config.json')
    edge_properties = [
        'afferent_section_id', 'afferent_segment_id', 'afferent_segment_offset',
        Synapse.U_SYN, Synapse.D_SYN, Synapse.F_SYN, Synapse.G_SYNX,
    ]
    morph_id = 110
    synapses = circuit.edges['default'].afferent_edges(morph_id, edge_properties)
    morph_filepath = circuit.nodes['All'].morph.get_filepath(morph_id)
    m = nm.load_morphology(morph_filepath)

    fig = morphology.draw(m, synapses)
    fig.show()


if __name__ == '__main__':
    # circuit_example()
    # example_afferent()
    example_efferent()
    # example_afferent_efferent()
