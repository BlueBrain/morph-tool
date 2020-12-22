"""Resampling for morphio morphologies
"""

import numpy as np
import morphio


def _vertex_path_lengths(points):
    """Using the segment length assign a path length for each vertex in the
    points array.

    Args:
        segment_lengths (np.ndarray): (N,) array of floats corresponding to
            segment lengths of a section

    Example:
        [1, 2, 3] -> [0, 1, 3, 6]

    Notes:
        The last element in the path length array corresponds to the total
        length of the segments.
    """
    segment_lengths = np.linalg.norm(points[1:] - points[:-1], axis=1)

    path_lengths = np.empty(len(segment_lengths) + 1, dtype=np.float32)
    path_lengths[0] = 0.
    path_lengths[1:] = np.cumsum(segment_lengths)
    return path_lengths


def _resample_from_linear_density(points, linear_density):
    """Creates a resampling of K-2 elements on the polyline defined by
    N `points`, omitting the first and last values which should remain
    unchanged. The number of elements is determined by the `linear_density`.

    The new sampling is represented by the line segment parametric equation:

    v'[k] = v[ids[k]] + fractions[k] * (v[ids[k] + 1] - v[ids[k]])

    where (v'[k], v[k]) can be either `points` or any other property mapped
    on `points`.

    Args:
        points (np.ndarray): (N, 3) Array of consecutive points definining
        segments.
        linear_density (float): Target number of points per micron

    Returns:
        tuple of numpy arrays:
            - ids (np.array): (K-2,) Array of positional ids
            - fractions (np.array): (K-2,) Parametric fractions

    Notes:
        Given consecutive points [p0, p1, ..., pN] a resampling will create
        new points on the polyline defined by these points. This means
        that a new point p' can be defined as a fraction t of an existing
        segment (p[n], p[n + 1]). Therefore, using the line segment's
        parametric eq:

        p' = p[n] + t * (p[n+1] - p[n])

        where p[n] is the start of the segment, (p[n+1] - p[n]) the segment's
        vector and t the fraction [0., 1.] of the vector to create the point.

        Therefore, this function returns K ids and fractions, where each id n
        and fraction t corespond to tehe equation for p' above. In vectorized
        form will then be:

        p'[k] = points[ids[k]] + fractions[k] * (p[ids[k] + 1] - p[ids[k]])

        For example if we have four points and the new point p'[k] lies
        between p[1] and p[2], and exactly at the middle:

        p[0] ---- p[1] --p'[k]-- p[2] ---- p[3]

        then the equation will become:

        p'[k] = p[1] + 0.5 * (p[2] - p[1])

        This interpolation relation holds for all values mapped onto points.
        For example the respective new diameters d' for i will be:

        d'[k] = d[ids[k]] + fractions[k] * (d[ids[k] + 1] - d[ids[k]])
    """
    path_lengths = _vertex_path_lengths(points)
    total_length = path_lengths[-1]

    # segment number according to linear density, minimum 1
    n_segments = max(1, int(total_length * linear_density))

    # split total length in n_segments equal parts of dl length
    dl = total_length / n_segments

    # vertex path lengths without including first and last
    new_path_lengths = np.arange(dl, total_length, dl)

    # ids of the starting points for the parametric equation
    ids = np.searchsorted(path_lengths, new_path_lengths) - 1

    # each fraction t corresponds to the fraction of the vector
    # in the parametric equation of the line segment p0 + t * (p1 - p0)
    fractions = (new_path_lengths - path_lengths[ids]) / (path_lengths[ids + 1] - path_lengths[ids])

    return ids, fractions


def _parametric_values(values, ids, fractions):
    """It creates a new array of values, where the start and end
    values remain unchanged and the in-between values are calculated
    from the parametric equation:

    v'[k] = v[ids[k]] + fractions[k] * (v[ids[k] + 1] - v[ids[k]])

    Args:
        values (np.ndarray): 1D or 2D array of values
        ids (np.ndarray): (N,) Array of ints
        fractions (np.ndarray): (N,) Array of floats

    Returns:
        new_values (np.ndarray): 1D or 2D array of length N + 2
    """
    n_new_values = len(ids) + 2

    if values.ndim > 1:
        fractions = fractions[:, np.newaxis]
        new_values = np.empty_like(
            values,
            shape=(n_new_values, values.shape[1])
        )
    else:
        new_values = np.empty_like(
            values,
            shape=(n_new_values,)
        )

    new_values[0] = values[0]
    new_values[-1] = values[-1]
    new_values[1:-1] = values[ids] + fractions * (values[ids + 1] - values[ids])

    return new_values


def _resample_section_neuron(section, linear_density):
    """Resample in-place the section data based on linear density.
    Args:
        section (Section): Mutable morphology's section
        linear_density (float): Linear density to determine the point number
    """
    points = section.points
    ids, fractions = _resample_from_linear_density(points, linear_density)
    section.points = _parametric_values(points, ids, fractions)
    section.diameters = _parametric_values(section.diameters, ids, fractions)


def _resample_section_astrocyte(section, linear_density):
    """Resample in-place the section data based on linear density.
    Args:
        section (Section): Mutable morphology's section
        linear_density (float): Linear density to determine the point number
    """
    points = section.points
    ids, fractions = _resample_from_linear_density(points, linear_density)
    section.points = _parametric_values(points, ids, fractions)
    section.diameters = _parametric_values(section.diameters, ids, fractions)
    section.perimeters = _parametric_values(section.perimeters, ids, fractions)


def _dispatch_section_function(cell_family):
    """Returns section function depending on the cell_family
    of the morphology
    """
    return {
        morphio.CellFamily.FAMILY_NEURON: _resample_section_neuron,
        morphio.CellFamily.FAMILY_GLIA: _resample_section_astrocyte
    }[cell_family]


def resample_linear_density(obj, linear_density):
    """Resample the number of points in morphology

    Args:
        obj: Morphology object, mutable or immutable
        linear_density (float): Number of points per micron

    Returns:
        morphio.mut.Morphology: Resampled morphology copy

    Notes:
        Resampling maintains cell topology, i.e. it
        does not modify branching and termination points

        A linear density of 1 corresponds to 1 point per um, 0.5 (1/2) to
        1 point per 2 um, 2.0 (2/1) to two points per 1 um, etc.
    """
    obj = morphio.mut.Morphology(obj)
    resample_function = _dispatch_section_function(obj.cell_family)

    for section in obj.iter():
        resample_function(section, linear_density)

    return obj
