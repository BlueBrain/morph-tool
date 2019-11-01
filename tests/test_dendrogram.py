import os

import neurom as nm
import numpy as np
import pandas as pd
from bluepy.v2.enums import Synapse
from nose.tools import assert_raises
from plotly.offline import plot

from morph_tool.dendrogram import draw

PATH = os.path.join(os.path.dirname(__file__), 'data')


def _get_neuron():
    filename = os.path.join(PATH, 'simple.swc')
    return nm.load_neuron(filename)


def _get_synapses(postGIDs):
    columns = [
        Synapse.PRE_GID, Synapse.POST_GID,
        Synapse.POST_SECTION_ID, Synapse.POST_SECTION_DISTANCE,
        Synapse.PRE_SECTION_ID, Synapse.PRE_SECTION_DISTANCE,
    ]
    synapses_number = 3
    data = [[0, postGID, 1, 10., 1, 10.]
            for _ in range(synapses_number) for postGID in postGIDs]
    synapse_ids = np.random.randint(100, size=synapses_number)
    index = pd.MultiIndex.from_product([postGIDs, synapse_ids])
    synapses = pd.DataFrame(data, index=index, columns=columns)
    return synapses


def test_explicit_postgid():
    postGID = 1
    fig = draw(_get_neuron(), _get_synapses([postGID, 2]), postGID)
    plot(fig, auto_open=False)


def test_implicit_valid_postgid():
    postGID = 1
    fig = draw(_get_neuron(), _get_synapses([postGID]))
    plot(fig, auto_open=False)


def test_implicit_invalid_postgid():
    postGID = 1
    with assert_raises(ValueError) as cm:
        draw(_get_neuron(), _get_synapses([postGID, 2]))
    assert 'neuron_gid' in cm.exception.args[0].lower()
