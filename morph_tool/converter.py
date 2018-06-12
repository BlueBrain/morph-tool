'''A morphology converter that tries to keep the soma surface equal'''
import os

from morphio import MorphologyVersion, SomaType
from morphio.mut import Morphology

import numpy as np
from numpy.linalg import eig, norm

np.set_printoptions(precision=19)


XYZ = slice(3)
X, Y, Z, R = 0, 1, 2, 3
cX, cY, cXY = np.s_[:, X], np.s_[:, Y], np.s_[:, :Z]


def contourcenter(xyz):
    '''python implementation of NEURON code: lib/hoc/import3d/import3d_sec.hoc
    '''
    POINTS = 101

    points = np.vstack((np.diff(xyz[:, [X, Y]], axis=0), (0, 0)))
    perim = np.cumsum(np.hstack(((0, ), norm(points, axis=1))))[:-1]

    d = np.linspace(0, perim[-1], POINTS)
    new_xyz = np.zeros((POINTS, 3))
    for i in range(3):
        new_xyz[:, i] = np.interp(x=d, xp=perim, fp=xyz[:, i])

    mean = np.mean(new_xyz, axis=0)

    return mean, new_xyz


def get_sides(points, major, minor):
    '''
    Circular permutation of the points so that the point with the largest
    coordinate along the major axis becomes the last point
    tobj = major.c.mul(d.x[i])  ###### uneeded? line 1191
    '''

    major_coord, minor_coord = np.dot(points, major), np.dot(points, minor)

    imax = np.argmax(major_coord)
    major_coord, minor_coord = np.roll(
        major_coord, -imax), np.roll(minor_coord, -imax)

    imin = np.argmin(major_coord)

    sides = [major_coord[:imin][::-1], major_coord[imin:]]
    rads = [minor_coord[:imin][::-1], minor_coord[imin:]]
    return sides, rads


def make_convex(sides, rads):
    '''Keep only points that make path convex'''
    def convex_idx(m):
        '''Return index to elements of m that make it convex

        Note: not efficient at the moment
        # now we have the two sides without the min and max points (rads[0]=0)
        # we hope both sides now monotonically increase, i.e. convex
        # make it convex

        '''
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
    '''this follows the function in
            lib/hoc/import3d/import3d_gui.hoc
       most of the comments are from there, so if you want to follow along, it should
       break up the function the same way
    '''
    # find the major axis of the ellipsoid that best fits the shape
    # assuming (falsely in general) that the center is the mean

    points = (points - mean)
    eigen_values, eigen_vectors = eig(np.dot(points.T, points))

    # To be consistent with NEURON eigen vector directions
    eigen_vectors *= -1

    idx = np.argmax(eigen_values)
    major = eigen_vectors[:, idx]
    # minor is normal and in xy plane
    idx = 3 - np.argmin(eigen_values) - np.argmax(eigen_values)
    minor = eigen_vectors[:, idx]
    minor[2] = 0

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


def from_swc(neuron, output_ext):
    '''Convert to SWC'''
    if output_ext == 'swc':
        return neuron

    if neuron.soma_type == SomaType.SOMA_CYLINDERS:
        direction = neuron.soma.points[-1] - neuron.soma.points[0]

        # 90 degree rotation along Z axis
        orthogonal = np.array([direction[1], -direction[0], 0])
        orthogonal /= np.linalg.norm(orthogonal)
        orthogonal = np.repeat(
            orthogonal[np.newaxis, :], len(neuron.soma.points), axis=0)
        contour_side1 = neuron.soma.points + \
            (orthogonal.T * neuron.soma.diameters / 2.).T
        contour_side2 = neuron.soma.points - \
            (orthogonal.T * neuron.soma.diameters / 2.).T
        contour_side2 = contour_side2[::-1]

        neuron.soma.points = np.vstack((contour_side1, contour_side2))
        neuron.soma.diameters = [0] * len(neuron.soma.points)

    elif not neuron.soma_type == SomaType.SOMA_NEUROMORPHO_THREE_POINT_CYLINDERS:
        if output_ext in ('asc', 'h5'):
            # We convert the cylinder into a disk of same surface
            # To preserve surface we need the disk radius must be twice the cylinder radius
            radius = neuron.soma.diameters[0]
            N = 20
            points = np.zeros((N, 3))
            phase = 2 * np.pi / (N - 1) * np.arange(N)
            points[:, 0] = radius * np.cos(phase)
            points[:, 1] = radius * np.sin(phase)
            points += neuron.soma.points[0]
            neuron.soma.points = points
            neuron.soma.diameters = np.repeat(radius, N)
    else:
        raise Exception(
            'A SWC morphology is not supposed to have a soma of type: {}'.format(
                neuron.soma_type))

    return neuron


def from_h5_or_asc(neuron, output_ext):
    '''Convert to ASC/H5'''
    if neuron.soma_type != SomaType.SOMA_SIMPLE_CONTOUR:
        raise Exception(
            'A H5 file morphology is not supposed to have a soma of type: {}'.format(
                neuron.soma_type))

    if output_ext == 'swc':
        mean, new_xyz = contourcenter(neuron.soma.points)
        neuron.soma.points, neuron.soma.diameters = contour2centroid(
            mean, new_xyz)
    return neuron


def run(input_file, outputfile):
    '''Run the appropriate converter'''
    neuron = Morphology(input_file)

    if neuron.soma_type == SomaType.SOMA_SINGLE_POINT:
        neuron.write(outputfile)
        return

    output_ext = os.path.splitext(outputfile)[1]
    if output_ext not in ('.swc', '.asc', '.h5'):
        raise Exception('Output file format should be one swc, asc or h5')
    output_ext = output_ext[1:]  # Remove the dot

    try:
        converter = {MorphologyVersion.MORPHOLOGY_VERSION_SWC_1: from_swc,
                     MorphologyVersion.MORPHOLOGY_VERSION_ASC_1: from_h5_or_asc,
                     MorphologyVersion.MORPHOLOGY_VERSION_H5_1: from_h5_or_asc,
                     MorphologyVersion.MORPHOLOGY_VERSION_H5_1_1: from_h5_or_asc,
                     MorphologyVersion.MORPHOLOGY_VERSION_H5_2: from_h5_or_asc
                     }[neuron.version]
    except KeyError:
        raise Exception(
            'No converter for morphology type: {}'.format(neuron.version))

    converter(neuron, output_ext).write(outputfile)
