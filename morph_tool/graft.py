'''What we call grafting is the process of taking one piece of a neuron and
putting it on another neuron.'''
import logging

import numpy as np
from morphio import SectionType, SomaType, Section as ImmutSection
from morphio.mut import Section


from morph_tool import NoAxonException, NoDendriteException, MorphToolException
from morph_tool.transform import translate


logger = logging.getLogger('morph_tool')


def find_axon(axon_or_donor_axon):
    """Find the root section which is an axon.

    If the number of axons found is greater than 1, then return the first axon.
    """

    if isinstance(axon_or_donor_axon, (ImmutSection, Section)):
        if axon_or_donor_axon.type == SectionType.axon:
            return axon_or_donor_axon
        else:
            raise MorphToolException('The section: {} is not an axon'.format(axon_or_donor_axon))

    axons = [root for root in axon_or_donor_axon.root_sections
             if root.type == SectionType.axon]

    if not axons:
        raise NoAxonException('No axon found!')

    if len(axons) > 1:
        logger.warning("Found multiple axons (%s), choosing the first one.", len(axons))

    return axons[0]


def _dendrites_mean_direction(neuron):
    dendrites = [neurite for neurite in neuron.root_sections if neurite.type != SectionType.axon]
    if not dendrites:
        raise NoDendriteException("Neuron has no dendrites")

    directions = np.array([dendrite.points[-1] - dendrite.points[0] for dendrite in dendrites])
    return np.mean(directions, axis=0)


def _axon_dendrites_angle(neuron):
    axon = find_axon(neuron)
    axon_dir = axon.points[-1] - axon.points[0]
    dendrite_mean_direction = _dendrites_mean_direction(neuron)
    cos_angle = dendrite_mean_direction.dot(
        axon_dir) / np.linalg.norm(dendrite_mean_direction) / np.linalg.norm(axon_dir)
    return np.arccos(cos_angle)


def _rotation_around_axis(axis, angle):
    """Returns a rotation matrix around the selected axis by an angle."""
    d = np.array(axis, dtype=np.float) / np.linalg.norm(axis)

    sn = np.sin(angle)
    cs = np.cos(angle)

    eye = np.eye(3, dtype=np.float)
    skew = np.array([[0, -d[2], d[1]],
                     [d[2], 0, -d[0]],
                     [-d[1], d[0], 0]], dtype=np.float)

    mtx = eye + sn * skew + (1. - cs) * np.linalg.matrix_power(skew, 2)
    return mtx


def _rotate_vector(vec, axis, angle):
    """Rotates the input vector vec
       by a selected angle
       around a specific axis.
    """
    return np.dot(_rotation_around_axis(axis, angle), vec)


def _random_direction(axis, angle):
    '''Returns an direction that makes an theta angle 'angle'
    with 'axis' and that has a random phi angle in [0, 2 pi]'''
    orthogonal = np.cross(axis, [0, 0, 1])
    if np.linalg.norm(orthogonal) < 1e-7:
        orthogonal = np.cross(axis, [0, 1, 0])

    orthogonal = _rotate_vector(orthogonal, axis, np.random.uniform(2 * np.pi))
    vec = _rotate_vector(axis, orthogonal, angle)
    vec /= np.linalg.norm(vec)
    return vec


def _section_initial_direction(section):
    '''Returns the section initial direction'''
    direction = section.points[1] - section.points[0]

    # In case of duplicate points, we skip first point
    if np.linalg.norm(direction) < 1e-8:
        return section.points[2] - section.points[1]
    return direction


def _soma_mean_radius(neuron, soma_center):
    '''Returns the soma mean radius'''
    if neuron.soma_type == SomaType.SOMA_SINGLE_POINT:
        return neuron.soma.diameters[0] / 2.
    else:
        return np.mean(np.linalg.norm(neuron.soma.points - soma_center, axis=1))


def grafting_position(neuron, axon_or_donor_neuron):
    '''Find out where to put the grafted axon

    If the neuron to be grafted already has a neuron reuse its starting point

    Otherwise try to mimic the position of the axon with respect to dendrites in
    the donor axon. The position is choosen in order to
    preserve the direction between the axon initial segment direction
    and all other dendrites initial segments directions

    Returns:
        The position where to graft the neuron
    '''
    try:
        old_axon = find_axon(neuron)
        return old_axon.points[0]
    except NoAxonException:
        axon_direction = _random_direction(axis=_dendrites_mean_direction(neuron),
                                           angle=_axon_dendrites_angle(axon_or_donor_neuron))
        soma_center = np.mean(neuron.soma.points, axis=0)
        axon_start = soma_center + _soma_mean_radius(neuron, soma_center) * axon_direction
        return axon_start


def delete_old_axons(neuron):
    '''Delete all axons'''
    for root_section in neuron.root_sections:
        if root_section.type == SectionType.axon:
            neuron.delete_section(root_section)


def graft_axon(neuron, axon_or_donor_neuron):
    '''Graft an axon
    args:
        neuron (morphio.mut.Morphology): Neuron where the axon will be grafted
        axon_or_donor_neuron (morphio.Morphology): An axon or a neuron with an axon

    This is a very simple implementation.
    No check are performed and the new axon is simply
    translated at the place where the old axon was.
    '''

    donor_axon = find_axon(axon_or_donor_neuron)
    position = grafting_position(neuron, axon_or_donor_neuron)

    delete_old_axons(neuron)

    new_axon = neuron.append_root_section(donor_axon, recursive=True)

    translate(new_axon, position - new_axon.points[0])
