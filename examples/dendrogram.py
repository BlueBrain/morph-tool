import numpy as np
import pandas as pd

import neurom as nm
from bluepysnap import Circuit
from bluepysnap.sonata_constants import Edge
from bluepysnap.bbp import Synapse

from morph_tool.dendrogram import draw


def plain_example():
    """Example that shows how to draw a neuron dendrogram with a plain synapses dataframe."""
    # Those properties are required in synapses dataframe for positioning
    required_synapse_properties = [
        Edge.SOURCE_NODE_ID, Edge.TARGET_NODE_ID,
        Edge.POST_SECTION_ID, Edge.POST_SECTION_POS,
        Edge.PRE_SECTION_ID, Edge.PRE_SECTION_POS,
    ]
    # if you don't want to use BluePySnap then you can use plain string names
    required_synapse_properties = [
        '@source_node', '@target_node',
        'afferent_section_id', 'afferent_section_pos',
        'efferent_section_id', 'efferent_section_pos',
    ]
    data = np.array([
        [0, 116, 158, 0.81408846, 3, 0.7344886],
        [0, 116, 101, 0.045983203, 294, 0.24454929],
        [0, 116, 128, 0.018469656, 295, 0.4290702],
        [116, 0, 380, 0.84480673, 1, 0.29180855],
        [116, 0, 380, 0.8815143, 1, 0.38261607],
        [116, 0, 380, 0.89492387, 1, 0.4126523],
        [116, 0, 380, 0.90750885, 1, 0.44434384],
        [116, 0, 380, 0.91967326, 1, 0.47111204],
        [116, 0, 380, 0.93142796, 1, 0.49867398],
    ])
    synapses = pd.DataFrame(columns=required_synapse_properties, data=data)
    # fix wrong Pandas dtypes autodetect.
    synapses = synapses.astype({'@target_node': int, '@source_node': int,
                                'afferent_section_id': int, 'efferent_section_id': int})
    # convert '@target_node' to string for discrete colormap
    synapses = synapses.astype({'@target_node': str})
    neuron = nm.load_neuron('dendrogram_plain_example.swc')
    fig = draw(neuron, synapses, 116)
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
        [0.29116565, 116.06434, 10.367496, 3.1585026],
        [0.29116565, 116.06434, 10.367496, 3.1585026],
        [0.29116565, 116.06434, 10.367496, 3.1585026],
        [0.29116565, 116.06434, 10.367496, 3.1585026],
    ])
    synapses_data = pd.DataFrame(columns=synapse_data_properties, data=data)
    synapses = pd.concat([synapses, synapses_data], axis=1)
    fig = draw(neuron, synapses, 116)
    fig.show()


def circuit_example():
    """Example that shows how to draw a neuron dendrogram with synapses from a bluepysnap circuit.

    To make this example work, you would need a proper SONATA circuit.
    """
    circuit = Circuit('/path/to/sonata_circuit_config.json')
    edge_properties = [
        # Must have properties names. They are required for positioning of synapses.
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
    neuron = nm.load_neuron(morph_filepath)

    fig = draw(neuron, synapses, target_node_id)
    fig.show()


if __name__ == '__main__':
    plain_example()
    # circuit_example()
