from pathlib import Path
import neurom as nm
import numpy as np
import pandas as pd
from bluepysnap.sonata_constants import Edge
from plotly.offline import plot

from nose.tools import assert_raises

from morph_tool import dendrogram

DATA = Path(__file__).resolve().parent / 'data'


def _create_test_neuron():
    return nm.load_neuron(DATA / 'simple.swc')


def _create_test_synapses(target_node_ids):
    columns = [
        Edge.SOURCE_NODE_ID, Edge.TARGET_NODE_ID,
        Edge.POST_SECTION_ID, Edge.POST_SECTION_POS,
        Edge.PRE_SECTION_ID, Edge.PRE_SECTION_POS,
    ]
    synapses_number = 5
    data = [[0, target_node_id, 1, 0.5, 1, 0.7]
            for _ in range(synapses_number) for target_node_id in target_node_ids]
    synapse_ids = np.arange(synapses_number)
    index = pd.MultiIndex.from_product([target_node_ids, synapse_ids])
    synapses = pd.DataFrame(data, index=index, columns=columns)
    return synapses


def test_constants():
    assert dendrogram.POST_SECTION_ID == Edge.POST_SECTION_ID
    assert dendrogram.POST_SECTION_POS == Edge.POST_SECTION_POS
    assert dendrogram.TARGET_NODE_ID == Edge.TARGET_NODE_ID
    assert dendrogram.PRE_SECTION_ID == Edge.PRE_SECTION_ID
    assert dendrogram.PRE_SECTION_POS == Edge.PRE_SECTION_POS
    assert dendrogram.SOURCE_NODE_ID == Edge.SOURCE_NODE_ID


def test_explicit_neuron_node_id():
    neuron_node_id = 1
    fig = dendrogram.draw(
        _create_test_neuron(), _create_test_synapses([neuron_node_id, 2]), neuron_node_id)
    # returns a string that contains the HTML <div>, without saving to file
    plot(fig, auto_open=False, output_type='div')


def test_implicit_valid_neuron_node_id():
    neuron_node_id = 1
    fig = dendrogram.draw(_create_test_neuron(), _create_test_synapses([neuron_node_id]))
    # returns a string that contains the HTML <div>, without saving to file
    plot(fig, auto_open=False, output_type='div')


def test_implicit_invalid_neuron_node_id():
    neuron_node_id = 1
    with assert_raises(ValueError) as cm:
        dendrogram.draw(_create_test_neuron(), _create_test_synapses([neuron_node_id, 2]))
    assert 'neuron_node_id' in cm.exception.args[0].lower()
