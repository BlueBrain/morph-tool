"""A morphology converter that tries to keep the soma surface equal."""
import logging
from pathlib import Path
import tempfile

import numpy as np
from scipy.optimize import minimize_scalar
from morphio import Option, SomaType
from morphio._morphio import WriterError  # pylint: disable=no-name-in-module
from morphio.mut import Morphology
from neurom import morphmath

from morph_tool import transform
from morph_tool.exceptions import MorphToolException

L = logging.getLogger(__name__)

XYZ = slice(3)
X, Y, Z, R = 0, 1, 2, 3
cX, cY, cXY = np.s_[:, X], np.s_[:, Y], np.s_[:, :Z]


def contourcenter(xyz):
    """Python implementation of NEURON code lib/hoc/import3d/import3d_sec.hoc."""
    POINTS = 101

    points = np.vstack((np.diff(xyz[:, [X, Y]], axis=0), xyz[0, [X, Y]]))
    perim = np.cumsum(np.hstack(((0, ), np.linalg.norm(points, axis=1))))[:-1]

    d = np.linspace(0, perim[-1], POINTS)
    new_xyz = np.zeros((POINTS, 3))
    for i in range(3):
        new_xyz[:, i] = np.interp(x=d, xp=perim, fp=xyz[:, i])

    mean = np.mean(new_xyz, axis=0)

    return mean, new_xyz


def get_sides(points, major, minor):
    """Circular permutation of the points.

    So that the point with the largest coordinate along the major axis becomes the last point
    ``tobj = major.c.mul(d.x[i])``  ###### unneeded? line 1191
    """
    major_coord, minor_coord = np.dot(points, major), np.dot(points, minor)

    imax = np.argmax(major_coord)
    # pylint: disable=invalid-unary-operand-type
    major_coord, minor_coord = (np.roll(major_coord, -imax),
                                np.roll(minor_coord, -imax))

    imin = np.argmin(major_coord)

    sides = [major_coord[:imin][::-1], major_coord[imin:]]
    rads = [minor_coord[:imin][::-1], minor_coord[imin:]]
    return sides, rads


def make_convex(sides, rads):
    """Keep only points that make path convex."""
    def convex_idx(m):
        """Return index to elements of m that make it convex.

        Note: not efficient at the moment
        # now we have the two sides without the min and max points (rads[0]=0)
        # we hope both sides now monotonically increase, i.e. convex
        # make it convex
        """
        idx = np.ones_like(m, dtype=bool)
        last_val = m[-1]
        for i in range(len(m) - 2, -1, -1):
            if m[i] < last_val:
                last_val = m[i]
            else:
                idx[i] = False
        return idx

    for i_side in [0, 1]:
        ci = convex_idx(sides[i_side])
        sides[i_side] = sides[i_side][ci]
        rads[i_side] = rads[i_side][ci]

    return sides, rads


def contour2centroid(mean, points):
    """This follows the function in lib/hoc/import3d/import3d_gui.hoc.

    Most of the comments are from there, so if you want to follow along, it should
    break up the function the same way.
    """
    L.info('Converting soma contour into a stack of cylinders')

    # find the major axis of the ellipsoid that best fits the shape
    # assuming (falsely in general) that the center is the mean

    points = points - mean
    eigen_values, eigen_vectors = np.linalg.eig(np.dot(points.T, points))

    # To be consistent with NEURON eigen vector directions
    eigen_vectors *= -1

    idx = np.argmax(eigen_values)
    major = eigen_vectors[:, idx]
    # minor is normal and in xy plane
    idx = 3 - np.argmin(eigen_values) - np.argmax(eigen_values)
    minor = eigen_vectors[:, idx]
    minor[2] = 0
    minor /= np.linalg.norm(minor)

    sides, rads = get_sides(points, major, minor)
    sides, rads = make_convex(sides, rads)

    tobj = np.sort(np.hstack(sides))
    new_major_coord = np.linspace(tobj[1], tobj[-2], 21)
    rads[0] = np.interp(new_major_coord, sides[0], rads[0])
    rads[1] = np.interp(new_major_coord, sides[1], rads[1])

    points = major * new_major_coord[:, np.newaxis] + mean
    diameters = np.abs(rads[0] - rads[1])

    # avoid 0 diameter ends
    diameters[0] = np.mean(diameters[:2])
    diameters[-1] = np.mean(diameters[-2:])

    return points, diameters


def _create_contour(radius, point_count=20, line_width=0.25):
    """Create a contour from an origin and radius.

    Note: The contour is created in the xy plan.

    within the corpus of ASC files, the mean is 0.21 for the diameter.
    this diameter is used for displaying the width of the contour line in
    Neurolucida, and 0.25 seems to be a nice display value
    """
    points = np.zeros((point_count, 3))
    phase = 2 * np.pi / point_count * np.arange(point_count)
    points[:, 0] = radius * np.sin(phase)
    points[:, 1] = radius * np.cos(phase)
    diameters = np.repeat(line_width, point_count)

    return points, diameters


def cylinder_to_cylindrical_contour(neuron):
    """We convert the cylinder into a circular contour that represents the same sphere."""
    L.info('Converting 3 point soma to spherical soma with same surface')
    assert neuron.soma_type == SomaType.SOMA_NEUROMORPHO_THREE_POINT_CYLINDERS
    radius = neuron.soma.diameters[0] / 2.
    points, diameters = _create_contour(radius)
    origin = neuron.soma.points[0]
    neuron.soma.points = points + origin
    neuron.soma.diameters = diameters
    neuron.soma.type = SomaType.SOMA_SIMPLE_CONTOUR


def cylinder_stack_to_cylindrical_contour(neuron):
    """Convert the cylinder into a simple contour that represents the same sphere."""
    L.info('Converting soma stack of cylinders into a contour in the XY plane')
    assert neuron.soma_type == SomaType.SOMA_CYLINDERS

    direction = neuron.soma.points[-1] - neuron.soma.points[0]
    i = 2
    while np.isclose(direction, 0).all():
        if i > len(neuron.soma.points):
            msg = "Can not convert SOMA_CYLINDERS if all points are equal"
            raise MorphToolException(msg)
        direction = neuron.soma.points[-i] - neuron.soma.points[0]
        i = i + 1

    # 90 degree rotation along Z axis
    orthogonal = np.array([direction[1], -direction[0], 0])
    norm = np.linalg.norm(orthogonal)
    if abs(norm) < 1e-10:
        # if we can't rotate along Z, we do along X
        orthogonal = np.array([0, direction[2], -direction[1]])
        norm = np.linalg.norm(orthogonal)

    orthogonal /= norm
    orthogonal = np.repeat(
        orthogonal[np.newaxis, :], len(neuron.soma.points), axis=0)
    contour_side1 = neuron.soma.points + (orthogonal.T * neuron.soma.diameters / 2.).T
    contour_side2 = neuron.soma.points - (orthogonal.T * neuron.soma.diameters / 2.).T
    contour_side2 = contour_side2[::-1]

    neuron.soma.points = np.vstack((contour_side1, contour_side2))
    neuron.soma.diameters = [0] * len(neuron.soma.points)
    neuron.soma.type = SomaType.SOMA_SIMPLE_CONTOUR


def single_point_sphere_to_circular_contour(neuron, ensure_NRN_area=False):
    """Single point soma to contour soma.

    Transform a single point soma that represents a sphere into a circular contour that represents
    the same sphere.

    With ensure_NRN_area=True, we perform a search on radius to ensure an equivalent
    soma surface area, as measured by NEURON with morph_tool.neuron_surface.get_NEURON_surface.
    """
    L.info('Converting 1-point soma (spherical soma) to circular contour '
           'representing the same sphere')
    assert neuron.soma_type == SomaType.SOMA_SINGLE_POINT

    swc_radius = neuron.soma.diameters[0] / 2.0

    def make_soma(radius):
        points, diameters = _create_contour(radius)
        origin = neuron.soma.points[0]
        neuron.soma.points = points + origin
        neuron.soma.diameters = diameters
        neuron.soma.type = SomaType.SOMA_SIMPLE_CONTOUR

    make_soma(swc_radius)

    if ensure_NRN_area:
        surf = 4.0 * np.pi * swc_radius**2

        # pylint: disable=import-outside-toplevel
        from morph_tool.neuron_surface import get_NEURON_surface

        with tempfile.NamedTemporaryFile(suffix='.asc') as fp:

            def cost(radius):
                make_soma(radius)
                neuron.write(fp.name)
                return abs(surf - get_NEURON_surface(fp.name))

            radius = minimize_scalar(
                cost, bounds=(0.8 * swc_radius, 1.2 * swc_radius), options={"xatol": 1e-2}
            ).x
            make_soma(radius)


def contour_to_cylinders(neuron):
    """Convert the simple contour into a cylinder."""
    mean, new_xyz = contourcenter(neuron.soma.points)
    neuron.soma.points, neuron.soma.diameters = contour2centroid(
        mean, new_xyz)
    neuron.soma.type = SomaType.SOMA_CYLINDERS


def soma_to_single_point(soma):
    """Surface preserving cylindrical soma to a single point sphere."""
    L.info('Converting soma to a single point sphere, while preserving the surface')
    neurom_points = np.hstack((soma.points, 0.5 * soma.diameters[:, None]))
    surface_area = sum(morphmath.segment_area(seg)
                       for seg in zip(neurom_points[1:], neurom_points[:-1]))

    soma.points = np.mean(soma.points, axis=0)[None, :]
    soma.diameters = [float((surface_area / np.pi) ** 0.5)]
    soma.type = SomaType.SOMA_SINGLE_POINT


def from_swc(neuron, output_ext, ensure_NRN_area=False):
    """Convert to SWC."""
    if output_ext == 'swc':
        return neuron

    if neuron.soma_type == SomaType.SOMA_CYLINDERS:
        cylinder_stack_to_cylindrical_contour(neuron)
    elif neuron.soma_type == SomaType.SOMA_NEUROMORPHO_THREE_POINT_CYLINDERS:
        cylinder_to_cylindrical_contour(neuron)
    elif neuron.soma_type == SomaType.SOMA_SINGLE_POINT:
        single_point_sphere_to_circular_contour(neuron, ensure_NRN_area=ensure_NRN_area)
    else:
        raise MorphToolException(
            f'A SWC morphology is not supposed to have a soma of type: {neuron.soma_type}')

    return neuron


def from_h5_or_asc(neuron, output_ext, ensure_NRN_area=False):
    """Convert from ASC/H5."""
    if neuron.soma_type == SomaType.SOMA_SINGLE_POINT:
        if output_ext == 'asc':
            single_point_sphere_to_circular_contour(neuron, ensure_NRN_area=ensure_NRN_area)
    elif neuron.soma_type == SomaType.SOMA_SIMPLE_CONTOUR:
        if output_ext == 'swc':
            contour_to_cylinders(neuron)

    if output_ext != "h5":
        for sec in neuron.iter():
            sec.perimeters = []

    return neuron


def choose_target_soma_type(neuron, output_ext):
    """Choose a target soma type based on the given output extension."""
    if output_ext in ["asc", "h5"]:
        target_type = SomaType.SOMA_SIMPLE_CONTOUR
    elif neuron.soma_type == SomaType.SOMA_SINGLE_POINT:
        target_type = SomaType.SOMA_SINGLE_POINT
    elif neuron.soma_type == SomaType.SOMA_SIMPLE_CONTOUR:
        target_type = SomaType.SOMA_CYLINDERS
    else:
        target_type = None
    return target_type


def convert_soma(neuron, output_ext, ensure_NRN_area=False):
    """Convert the soma of the given neuron to the required soma type."""
    current_type = neuron.soma_type
    target_type = choose_target_soma_type(neuron, output_ext)
    L.info('Original soma type: %s', current_type)

    if target_type is not None and current_type != target_type:
        if target_type == SomaType.SOMA_SINGLE_POINT:
            soma_to_single_point(neuron.soma)
        elif target_type == SomaType.SOMA_CYLINDERS:
            contour_to_cylinders(neuron)
        elif target_type == SomaType.SOMA_SIMPLE_CONTOUR:
            if current_type == SomaType.SOMA_SINGLE_POINT:
                single_point_sphere_to_circular_contour(neuron, ensure_NRN_area=ensure_NRN_area)
            elif current_type == SomaType.SOMA_CYLINDERS:
                cylinder_stack_to_cylindrical_contour(neuron)
            elif current_type == SomaType.SOMA_NEUROMORPHO_THREE_POINT_CYLINDERS:
                cylinder_to_cylindrical_contour(neuron)


def convert(input_file,
            output_file,
            recenter=False,
            nrn_order=False,
            single_point_soma=False,
            sanitize=False,
            ensure_NRN_area=False):
    """Run the appropriate converter.

    Args:
        input_file(str): path to input file
        output_file(str): path to output file
        recenter(bool): whether to recenter the morphology based on the
        center of gravity of the soma
        nrn_order(bool): whether to traverse the neuron in the NEURON fashion
        single_point_soma(bool): For SWC only
        sanitize(bool): whether to sanitize the morphology
        ensure_NRN_area(bool): to ensure area is preserved in NEURON from swc point soma
    """
    kwargs = {}
    if nrn_order:
        kwargs['options'] = Option.nrn_order

    neuron = Morphology(input_file, **kwargs)
    if sanitize:
        neuron.remove_unifurcations()

    output_ext = Path(output_file).suffix.lower()

    if output_ext not in ('.swc', '.asc', '.h5', ):
        raise MorphToolException('Output file format should be one swc, asc or h5')

    output_ext = output_ext[1:]  # Remove the dot

    convert_soma(neuron, output_ext, ensure_NRN_area=ensure_NRN_area)

    if single_point_soma:
        soma_to_single_point(neuron.soma)

    if output_ext != "h5":
        for sec in neuron.iter():
            sec.perimeters = []

    if recenter:
        transform.translate(neuron, -1 * neuron.soma.center)

    try:
        neuron.write(output_file)
    except WriterError as e:
        raise MorphToolException('Use `sanitize` option for converting') from e

    try:
        # pylint: disable=import-outside-toplevel
        from morph_tool.neuron_surface import get_NEURON_surface
        L.info('Soma surface as computed by NEURON:\n'
               'before conversion: %s\n'
               'after conversion: %s',
               get_NEURON_surface(input_file),
               get_NEURON_surface(output_file))
    except:  # noqa pylint: disable=bare-except
        L.info('Final NEURON soma surface check was skipped probably because BluePyOpt'
               ' or NEURON is not installed')

    return neuron
