'''Dendrogram helper functions and class'''
import numpy as np
import plotly.express as px
from plotly import graph_objects
from neurom import NeuriteType
from neurom.core import Neurite, Neuron
from neurom.view.dendrogram import Dendrogram, layout_dendrogram, move_positions, get_size
from neurom.view.view import TREE_COLOR

POST_SECTION_ID = 'afferent_section_id'
POST_SECTION_POS = 'afferent_section_pos'
TARGET_NODE_ID = '@target_node'
PRE_SECTION_ID = 'efferent_section_id'
PRE_SECTION_POS = 'efferent_section_pos'
SOURCE_NODE_ID = '@source_node'


class SynDendrogram(Dendrogram):
    """Dendrogram that keeps track of ``neurom_section.id`` that was used to create it."""

    def __init__(self, neurom_section):
        super().__init__(neurom_section)
        if isinstance(neurom_section, Neuron):
            self.section_id = -1
            self.children = [
                SynDendrogram(neurite.root_node) for neurite in neurom_section.neurites]
        else:
            if isinstance(neurom_section, Neurite):
                neurom_section = neurom_section.root_node
            self.section_id = neurom_section.id
            self.children = [SynDendrogram(child) for child in neurom_section.children]


def _draw_polygon(coords, color):
    path = 'M' + ' L'.join('{:.2f},{:.2f}'.format(coord[0], coord[1]) for coord in coords) + ' Z'

    return graph_objects.layout.Shape(
        opacity=0.4, type="path", path=path, fillcolor=color, line_color=color)


def _draw_line(start, end, color):
    return graph_objects.layout.Shape(
        type="line", xref="x", yref="y", x0=start[0], y0=start[1], x1=end[0], y1=end[1],
        opacity=0.4,
        line=dict(
            color=color,
            width=2,
        ),
    )


def _get_draw_shapes(dendrogram, positions):
    """Generates drawable shapes for dendrogram.

    Args:
        dendrogram (Dendrogram): dendrogram
        positions (dict of Dendrogram: np.array): xy coordinates of dendrograms

    Returns:
        list: plotly shapes
    """
    color = TREE_COLOR[dendrogram.neurite_type]
    shapes = [_draw_polygon(dendrogram.coords + positions[dendrogram], color)]
    if dendrogram.children:
        end_point = positions[dendrogram] + [0, dendrogram.height]
        for child in dendrogram.children:
            shapes.append(_draw_line(end_point, positions[child], color))
            shapes += _get_draw_shapes(child, positions)
    return shapes


def _position_synapses(positions, synapses, neuron_node_id):
    """Position synapses according to positions of dendrogram.

    Adjusts ``synapses`` inplace with new columns for xy coordinates.

    Args:
        positions (dict of SynDendrogram: np.array): positions xy coordinates of dendrograms
        synapses (pd.DataFrame): synapses dataframe
        neuron_node_id: neuron's node id from which dendrogram has been built
    """

    def _jitter_x(grp):
        if len(grp) > 1:
            min_x = grp['x'].min()
            jittered_x = np.arange(0, len(grp), 1) + min_x
            grp.loc[:, 'x'] = jittered_x
        return grp

    for dendrogram, position in positions.items():
        post_section = ((synapses[POST_SECTION_ID] == dendrogram.section_id)
                        & (synapses[TARGET_NODE_ID] == neuron_node_id))
        synapses.loc[post_section, 'x'] = position[0]
        synapses.loc[post_section, 'y'] = position[1] + synapses.loc[
            post_section, POST_SECTION_POS] * dendrogram.height
        pre_section = ((synapses[PRE_SECTION_ID] == dendrogram.section_id)
                       & (synapses[SOURCE_NODE_ID] == neuron_node_id))
        synapses.loc[pre_section, 'x'] = position[0]
        synapses.loc[pre_section, 'y'] = position[1] + synapses.loc[
            pre_section, PRE_SECTION_POS] * dendrogram.height

        # jitter too close synapses
        df = synapses[pre_section | post_section]
        if len(df) > 1:
            # pylint: disable=cell-var-from-loop
            # group by Y coordinate and jitter by X coordinate
            synapses[pre_section | post_section] = df \
                .groupby(lambda idx: round(df.loc[idx, 'y'])) \
                .apply(_jitter_x)


def _add_neurite_legend(fig, dendrogram):
    def neurite_legend(neurite_type):
        fig.add_trace(graph_objects.Scatter(
            x=[None], y=[None],
            name=neurite_type.name,
            line={'color': TREE_COLOR[neurite_type]}
        ))

    if dendrogram.neurite_type != NeuriteType.soma:
        neurite_legend(dendrogram.neurite_type)
    else:
        neurite_types = {d.neurite_type for d in [dendrogram] + dendrogram.children}
        for neurite_type in neurite_types:
            neurite_legend(neurite_type)


def _get_default_neuron_node_id(synapses):
    node_id_set = synapses[TARGET_NODE_ID].unique()
    if node_id_set.size > 1:
        raise ValueError('Please specify `neuron_node_id` parameter.')
    return node_id_set[0]


def _get_synapse_colormap(synapses):
    color_list = [
        'rgb(0,100,0)', 'rgb(139,69,19)', 'rgb(112,128,144)', 'rgb(255,215,0)', 'rgb(218,112,214)']
    node_id_set = synapses[TARGET_NODE_ID].unique()
    return {id_: color_list[idx % len(color_list)] for idx, id_ in enumerate(node_id_set)}


def draw(neuron, synapses=None, neuron_node_id=None):
    """Draw dendrogram with synapses.

    Args:
        neuron (Neurite|Neuron): a Neurite or a Neuron instance of NeuroM package.
        synapses (DataFrame): synapses dataframe.
        neuron_node_id (int|None): node id of ``neuron``. If None then it is taken from
            ``synapses[TARGET_NODE_ID]``.

    Returns:
        plotly.graph_objects.Figure: plotly figure
    """
    synapses = synapses.astype({'@target_node': int, '@source_node': int,
                                'afferent_section_id': int, 'efferent_section_id': int})
    dendrogram = SynDendrogram(neuron)
    positions = layout_dendrogram(dendrogram, np.array([0, 0]))
    w, h = get_size(positions)
    positions = move_positions(positions, np.array([.5 * w, 0]))
    if not synapses.empty:
        if neuron_node_id is None:
            neuron_node_id = _get_default_neuron_node_id(synapses)
        _position_synapses(positions, synapses, neuron_node_id)
        synapses['synapse_id'] = synapses.index.values
        # convert 'target node' to string for discrete colormap
        synapses = synapses.astype({TARGET_NODE_ID: str})
        fig = px.scatter(
            synapses, x='x', y='y', color=TARGET_NODE_ID,
            color_discrete_map=_get_synapse_colormap(synapses),
            hover_data=synapses.columns)
    else:
        fig = graph_objects.Figure()
    fig.update_xaxes(range=[-.05 * w, 1.05 * w], zeroline=False)
    fig.update_yaxes(range=[-.05 * h, 1.05 * h], zeroline=False)
    shapes = _get_draw_shapes(dendrogram, positions)
    fig.update_layout(shapes=shapes)
    _add_neurite_legend(fig, dendrogram)
    return fig
