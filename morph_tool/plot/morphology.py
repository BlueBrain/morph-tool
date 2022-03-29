"""Module for drawing morphologies with synapses."""

from functools import partial
import warnings

import pandas as pd
from pandas import Categorical
from pandas.core.dtypes.common import is_integer_dtype
from plotly import express

from neurom import morphmath, COLS
from neurom.view.plotly_impl import get_figure
from morph_tool.plot.consts import (PRE_SECTION_ID, PRE_SEGMENT_ID, PRE_SEGMENT_OFFSET,
                                    POST_SECTION_ID, POST_SEGMENT_ID, POST_SEGMENT_OFFSET)

REQUIRED_PRE_COLUMNS = {PRE_SECTION_ID, PRE_SEGMENT_ID, PRE_SEGMENT_OFFSET}
REQUIRED_POST_COLUMNS = {POST_SECTION_ID, POST_SEGMENT_ID, POST_SEGMENT_OFFSET}


def _validate_synapses(synapses):
    def _is_int(*columns):
        return all(is_integer_dtype(synapses[c].dtype) for c in columns)

    is_pre = REQUIRED_PRE_COLUMNS.issubset(synapses.columns)
    is_post = REQUIRED_POST_COLUMNS.issubset(synapses.columns)
    assert is_pre or is_post, \
        f'{REQUIRED_PRE_COLUMNS} or {REQUIRED_POST_COLUMNS} are required columns of `synapses`'
    if is_pre and not _is_int(PRE_SECTION_ID, PRE_SEGMENT_ID) or \
            is_post and not _is_int(POST_SECTION_ID, POST_SEGMENT_ID):
        warnings.warn('Section ids and segment ids columns of `synapses` arg are not integers, and '
                      'will be forced to integer type')


def _add_coords(synapse, morphology):
    """Adds coordinates and direction fields to ``synapses`` via ``apply`` function."""
    is_pre = PRE_SECTION_ID in synapse.index and not pd.isnull(synapse[PRE_SECTION_ID])
    is_post = POST_SECTION_ID in synapse.index and not pd.isnull(synapse[POST_SECTION_ID])
    assert is_pre != is_post, "Synapse must have either afferent or efferent section ids for the " \
                              "morphology. It can't have both at the same time."

    if is_post:
        sec_id = int(synapse[POST_SECTION_ID])
        seg_id = int(synapse[POST_SEGMENT_ID])
        seg_ofst = float(synapse[POST_SEGMENT_OFFSET])
        synapse['direction'] = 'afferent'
    else:
        sec_id = int(synapse[PRE_SECTION_ID])
        seg_id = int(synapse[PRE_SEGMENT_ID])
        seg_ofst = float(synapse[PRE_SEGMENT_OFFSET])
        synapse['direction'] = 'efferent'

    if sec_id == 0:
        # synapse is on soma
        p = morphology.soma.points[0]
        synapse['x'], synapse['y'], synapse['z'] = p[COLS.XYZ]
        # place synapse on surface of soma along Z axes so it won't be hidden by soma on the plot
        synapse['z'] += p[COLS.R]
        return synapse

    # NeuroM morphology counts sections from 0. Synapses count sections from 1 because section
    # id 0 is for soma.
    sec_id -= 1
    sec = morphology.sections[sec_id]
    assert sec_id == sec.id, f'Error. Synapse with section id {sec_id} must map to the same ' \
                             f'section id in `morphology` arg but maps to {sec.id}.'
    assert 0 <= seg_id <= len(sec.points), f'No such segment id {seg_id} for section id ' \
                                           f'{sec_id} of `morphology` arg'

    seg_p1, seg_p2 = sec.points[seg_id - 1], sec.points[seg_id]
    seg_len = morphmath.point_dist(seg_p1, seg_p2)
    coords = morphmath.linear_interpolate(seg_p1, seg_p2, seg_ofst / seg_len)
    synapse['x'], synapse['y'], synapse['z'] = coords
    return synapse


def draw(morphology, synapses):
    """Draw morphology with synapses.

    Args:
        morphology (Neuron): a Neuron instance of NeuroM package.
        synapses (DataFrame): synapses dataframe. It is required to have columns from
          :py:mod:`morph_tool.plot.consts`. See an example for details.

    Returns:
        plotly.graph_objects.Figure: plotly figure

    Afferent only synapses example
        .. literalinclude:: ../../../examples/draw_morphology.py
            :pyobject: example_afferent

    Afferent and efferent synapses example
        .. literalinclude:: ../../../examples/draw_morphology.py
            :pyobject: example_afferent_efferent

    Circuit synapses example
        .. literalinclude:: ../../../examples/draw_morphology.py
            :pyobject: circuit_example

    """
    _validate_synapses(synapses)
    synapses['direction'] = Categorical(['None'] * len(synapses),
                                        categories=['afferent', 'efferent', 'None'])
    synapses = synapses.apply(partial(_add_coords, morphology=morphology), axis=1)

    fig = express.scatter_3d(synapses, x='x', y='y', z='z',
                             color='direction',
                             color_discrete_map={'afferent': 'orange', 'efferent': 'green',
                                                 'None': 'black'},
                             hover_data=synapses.columns)

    figure_dict = get_figure(morphology, plane='3d', title=morphology.name)
    fig.add_traces(figure_dict['data'])
    fig.update_layout(figure_dict['layout'])
    return fig
