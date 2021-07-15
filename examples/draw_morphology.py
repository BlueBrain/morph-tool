import numpy as np
import pandas as pd
import neurom as nm
from morph_tool.plot import morphology
from morph_tool.plot import consts


def example_afferent():
    neuron = nm.load_neuron('dendrogram_plain_example.swc')
    data = np.array([
        [4, 0, 0.81408846, 'additional value'],
        [5, 0, 0.145983203, 'additional value'],
        [1, 0, 0.4290702, 'additional value'],
    ])
    columns = [consts.POST_SECTION_ID, consts.POST_SEGMENT_ID, consts.POST_SEGMENT_OFFSET,
               'additional_data']
    synapses = pd.DataFrame(data, columns=columns)

    fig = morphology.draw(neuron, synapses)
    fig.show()


def example_efferent():
    neuron = nm.load_neuron('dendrogram_plain_example.swc')
    data = np.array([
        [4, 0, 0.81408846, 'additional value'],
        [5, 0, 0.145983203, 'additional value'],
        [1, 0, 0.4290702, 'additional value'],
    ])
    columns = [consts.PRE_SECTION_ID, consts.PRE_SEGMENT_ID, consts.PRE_SEGMENT_OFFSET,
               'additional_data']
    synapses = pd.DataFrame(data, columns=columns)

    fig = morphology.draw(neuron, synapses)
    fig.show()


def example_afferent_efferent():
    neuron = nm.load_neuron('dendrogram_plain_example.swc')
    data = np.array([
        [4, 0, 0.81408846, np.nan, np.nan, np.nan],
        [5, 0, 0.145983203, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, 1, 0, 0.4290702],
        [np.nan, np.nan, np.nan, 1, 0, 0.29180855],
        [np.nan, np.nan, np.nan, 1, 0, 0.68261607],
    ])
    columns = [consts.POST_SECTION_ID, consts.POST_SEGMENT_ID, consts.POST_SEGMENT_OFFSET,
               consts.PRE_SECTION_ID, consts.PRE_SEGMENT_ID, consts.PRE_SEGMENT_OFFSET]
    synapses = pd.DataFrame(data, columns=columns)

    fig = morphology.draw(neuron, synapses)
    fig.show()


if __name__ == '__main__':
    example_afferent()
    # example_efferent()
    # example_afferent_efferent()
