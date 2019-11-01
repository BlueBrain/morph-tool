'''Dendrogram helper functions and class'''
import numpy as np
import plotly.express as px
from plotly import graph_objects
from bluepy.v2 import Synapse
from neurom import NeuriteType
from neurom.core import Neurite, Neuron
from neurom.view.dendrogram import Dendrogram, layout_dendrogram, move_positions, get_size
from neurom.view.view import TREE_COLOR


class SynDendrogram(Dendrogram):
    '''SynDendrogram. Dendrogram that keeps track of neurom_section.id that was used to create it.
    '''
    def __init__(self, neurom_section):
        super(SynDendrogram, self).__init__(neurom_section)
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
    path = 'M' + ' L'.join(['{:.2f},{:.2f}'.format(coord[0], coord[1]) for coord in coords]) + ' Z'

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
    '''Generates drawable shapes for dendrogram.

    Args:
        dendrogram (Dendrogram): dendrogram
        positions (dict of Dendrogram: np.array): positions xy coordinates of dendrograms

    Returns:
        List of plotly Shapes.
    '''
    color = TREE_COLOR[dendrogram.neurite_type]
    shapes = [_draw_polygon(dendrogram.coords + positions[dendrogram], color)]
    if dendrogram.children:
        end_point = positions[dendrogram] + [0, dendrogram.height]
        for child in dendrogram.children:
            shapes.append(_draw_line(end_point, positions[child], color))
            shapes += _get_draw_shapes(child, positions)
    return shapes


def _position_synapses(positions, synapses, neuron_gid):
    '''Position synapses according to positions of dendrogram.

    Adjusts ``synapses`` inplace with new columns for xy coordinates.
    Args:
        positions (dict of SynDendrogram: np.array): positions xy coordinates of dendrograms
        synapses: synapses dataframe from bluepy
        neuron_gid: neuron's gid from which dendrogram has been built
    '''

    def _jitter_x(grp):
        if len(grp) > 1:
            min_x = grp['x'].min()
            jittered_x = np.arange(0, len(grp), 1) + min_x
            grp.loc[:, 'x'] = jittered_x
        return grp

    for dendrogram, position in positions.items():
        post_section = ((synapses[Synapse.POST_SECTION_ID.value] == dendrogram.section_id)
                        & (synapses[Synapse.POST_GID.value] == neuron_gid))
        synapses.loc[post_section, 'x'] = position[0]
        synapses.loc[post_section, 'y'] = position[1] + synapses.loc[
            post_section, Synapse.POST_SECTION_DISTANCE.value]
        pre_section = ((synapses[Synapse.PRE_SECTION_ID.value] == dendrogram.section_id)
                       & (synapses[Synapse.PRE_GID.value] == neuron_gid))
        synapses.loc[pre_section, 'x'] = position[0]
        synapses.loc[pre_section, 'y'] = position[1] + synapses.loc[
            pre_section, Synapse.PRE_SECTION_DISTANCE.value]

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


def _get_default_neuron_gid(synapses):
    post_gid_set = synapses[Synapse.POST_GID.value].unique()
    if post_gid_set.size > 1:
        raise ValueError('Please specify `neuron_gid` parameter.')
    return post_gid_set[0]


def _get_synapse_colormap(synapses):
    color_list = [
        'rgb(0,100,0)', 'rgb(139,69,19)', 'rgb(112,128,144)', 'rgb(255,215,0)', 'rgb(218,112,214)']
    gid_list = synapses[Synapse.POST_GID.value].unique()
    return {gid: color_list[i % len(color_list)] for i, gid in enumerate(gid_list)}


def draw(neuron, synapses=None, neuron_gid=None):
    ''' Draw dendrogram with synapses.

    Args:
        neuron (Neurite|Neuron): a Neurite or a Neuron instance of NeuroM package.
        synapses (DataFrame): synapses dataframe from bluepy.
        neuron_gid (int|None): ``neuron`` gid. If None then it is taken from POST_GID.

    Returns:
        Plotly Figure.
    '''
    dendrogram = SynDendrogram(neuron)
    positions = layout_dendrogram(dendrogram, np.array([0, 0]))
    w, h = get_size(positions)
    positions = move_positions(positions, np.array([.5 * w, 0]))
    if not synapses.empty:
        # Rename columns cause plotly does not accept enums as names for DataFrame
        synapses = synapses.rename(lambda x: x.value, axis=1)
        if neuron_gid is None:
            neuron_gid = _get_default_neuron_gid(synapses)
        _position_synapses(positions, synapses, neuron_gid)
        synapses['synapse_id'] = synapses.index.get_level_values(1)
        # We can't use synapses[Synapse.POST_GID.value] as is for colors because it is numeric
        # and plotly uses continuous color range for numeric values.
        synapses[Synapse.POST_GID.value] = synapses[Synapse.POST_GID.value].astype(str)
        fig = px.scatter(
            synapses, x='x', y='y', color=Synapse.POST_GID.value,
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
