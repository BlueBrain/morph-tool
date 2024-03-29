import numpy as np
import pandas as pd
import neurom as nm
from morph_tool.plot import dendrogram, consts


def plain_example():
    """Example that shows how to plot a neuron dendrogram with a plain synapses dataframe."""
    # Those properties are required in synapses dataframe for positioning
    required_synapse_properties = [
        consts.SOURCE_NODE_ID, consts.TARGET_NODE_ID,
        consts.POST_SECTION_ID, consts.POST_SECTION_POS,
        consts.PRE_SECTION_ID, consts.PRE_SECTION_POS,
    ]
    # or use plain strings
    # required_columns = ['@source_node', '@target_node',
    #                     'afferent_section_id', 'afferent_section_pos',
    #                     'efferent_section_id', 'efferent_section_pos']
    data = np.array([
        [0, 116, 4, 0.81408846, 3, 0.7344886],
        [0, 116, 5, 0.145983203, 4, 0.24454929],
        [0, 116, 3, 0.968469656, 1, 0.4290702],
        [116, 0, 2, 0.84480673, 1, 0.29180855],
        [116, 0, 2, 0.5815143, 1, 0.68261607],
    ])
    synapses = pd.DataFrame(columns=required_synapse_properties, data=data)
    synapses = synapses.astype({'@target_node': int, '@source_node': int,
                                'afferent_section_id': int, 'efferent_section_id': int})
    m = nm.load_morphology('dendrogram_plain_example.swc')
    fig = dendrogram.draw(m, synapses, 116)
    fig.show()

    # If you want to show additional data with synapses then just use additional columns in you
    # dataframe. These data properties can have any names.
    synapse_data_properties = [
        'u_syn', 'depression_time', 'facilitation_time', 'conductance'
    ]
    data = np.array([
        [0.16718547, 153.8097, 8.452671, 1.9140357],
        [0.16718547, 153.8097, 8.452671, 1.9140357],
        [0.16718547, 153.8097, 8.452671, 1.9140357],
        [0.29116565, 116.06434, 10.367496, 3.1585026],
        [0.29116565, 116.06434, 10.367496, 3.1585026],
    ])
    synapses_data = pd.DataFrame(columns=synapse_data_properties, data=data)
    synapses = pd.concat([synapses, synapses_data], axis=1)
    fig = dendrogram.draw(m, synapses, 116)
    fig.show()


def circuit_example():
    """Example that shows how to plot a neuron dendrogram with synapses from a bluepysnap circuit.

    To make this example work, you would need a proper SONATA circuit.
    """
    from bluepysnap import Circuit
    from bluepysnap.sonata_constants import Edge
    from bluepysnap.bbp import Synapse
    circuit = Circuit('/path/to/sonata_circuit_config.json')
    # you can use `bluepysnap.sonata_constants.Edge` instead of `morph_tool.draw.consts`
    edge_properties = [
        Edge.SOURCE_NODE_ID, Edge.TARGET_NODE_ID,
        Edge.POST_SECTION_ID, Edge.POST_SECTION_POS,
        Edge.PRE_SECTION_ID, Edge.PRE_SECTION_POS,
        Synapse.U_SYN, Synapse.D_SYN, Synapse.F_SYN, Synapse.G_SYNX,
    ]
    source_node_id = 0
    target_node_id = 116
    synapses1 = circuit.edges['default'].pair_edges(source_node_id, target_node_id, edge_properties)
    synapses2 = circuit.edges['default'].pair_edges(target_node_id, source_node_id, edge_properties)
    synapses = pd.concat([synapses1, synapses2])

    morph_filepath = circuit.nodes['default'].morph.get_filepath(target_node_id)
    m = nm.load_morphology(morph_filepath)

    fig = dendrogram.draw(m, synapses, target_node_id)
    fig.show()


if __name__ == '__main__':
    plain_example()
    # circuit_example()
