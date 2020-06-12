from pathlib import Path
import neurom as nm
import numpy as np
import pandas as pd
from bluepysnap.sonata_constants import Edge
from nose.tools import assert_raises
from plotly.offline import plot

from morph_tool.dendrogram import draw

DATA = Path(__file__).resolve().parent / 'data'


def _create_test_neuron():
    return nm.load_neuron(DATA / 'simple.swc')


def _create_test_synapses(target_node_ids):
    columns = [
        Edge.SOURCE_NODE_ID, Edge.TARGET_NODE_ID,
        Edge.POST_SECTION_ID, Edge.POST_SECTION_POS,
        Edge.PRE_SECTION_ID, Edge.PRE_SECTION_POS,
    ]
    synapses_number = 3
    data = [[0, target_node_id, 1, 0.5, 1, 0.7]
            for _ in range(synapses_number) for target_node_id in target_node_ids]
    synapse_ids = np.random.randint(100, size=synapses_number)
    index = pd.MultiIndex.from_product([target_node_ids, synapse_ids])
    synapses = pd.DataFrame(data, index=index, columns=columns)
    return synapses


def test_explicit_neuron_node_id():
    neuron_node_id = 1
    fig = draw(_create_test_neuron(), _create_test_synapses([neuron_node_id, 2]), neuron_node_id)
    plot(fig, auto_open=False)


def test_implicit_valid_neuron_node_id():
    neuron_node_id = 1
    fig = draw(_create_test_neuron(), _create_test_synapses([neuron_node_id]))
    plot(fig, auto_open=False)


def test_implicit_invalid_neuron_node_id():
    neuron_node_id = 1
    with assert_raises(ValueError) as cm:
        draw(_create_test_neuron(), _create_test_synapses([neuron_node_id, 2]))
    assert 'neuron_node_id' in cm.exception.args[0].lower()
