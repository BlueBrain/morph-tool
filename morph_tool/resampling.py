"""Resampling functionality for morphio morphologies. It allows for the
generation of new points on the morphology skeleton with different density
"""
import numpy as np
import morphio


def _accumulated_path_lengths(points):
    """Return the accumulated path lengths along the polygonal line joining
    the consecutive points, starting with the point of index 0

    Args:
        points (np.ndarray): (N, 3) array of 3D points

    Returns:
        path_lengths (np.ndarray): (N,) Accumulated path lengths

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
    """Given a polyline of consecutive points that form segments, it creates
    a new sampling of K points (or point properties) on the same polyline. K
    is determined from the linear density as points per micron.

    Each new point / property is determined via two values: a starting point
    at positional id i in the points array i and a fraction along the segment formed
    by the consecutive points (points[i], points[i+1]). Thus, we have:

    v'[k] = v[ids[k]] + fractions[k] * (v[ids[k] + 1] - v[ids[k]])

    where ids are the positional indices in the points array, corresponding to the
    starting points to which fractions of the segments (v[ids[k]], v[ids[k+1]]) are
    added to reconstruct the new point or property.

    The reason of returning the ids and fractions instead of new points is because
    if instead of points a property is used, such as diameter, it will result to its
    iterpolation on the polyline. Therefore, the result of this function can be used
    to interpolate both points and point properties on a polyline.

    Args:
        points (np.ndarray): (N, 3) Array of consecutive points defining
        segments.
        linear_density (float): Target number of points per micron

    Returns:
        tuple:
            - ids (np.array): (K-2,) int array of positional ids
            - fractions (np.array): (K-2,) float array of fractions in the
                range [0, 1]

        Using the ids and fractions, the new interior points (without first
        and last) can be created:

        new_points = points[ids] + fractions * (points[ids + 1] - points[ids])

        Similarly, point properties (e.g. diameters) can be interpolated:

        new_diams = diams[ids] + fractions * (diams[ids + 1] - diams[ids])

    Notes:

        It doesn't not return the ids and fractions for first and last points in
        the new resampling, because they should remain unchanged. Therefore, the
        total number of ids and fractions is K-2.
    """
    path_lengths = _accumulated_path_lengths(points)
    total_length = path_lengths[-1]

    # segment number according to linear density, minimum 1
    # the number of points K = n_segments + 1
    n_segments = max(1, int(total_length * linear_density))

    # split total length in n_segments of equal length dl
    dl = total_length / n_segments

    # vertex path lengths without including first and last (K-2)
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

    See _resample_from_linear_density for more details on ids and fractions

    Args:
        values (np.ndarray): 1D or 2D array of values
            Points or any property defined on points can be interpolated
            using this function.
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


def _resample_neuron_section(section, linear_density):
    """Resample in-place the section data based on linear density.
    Args:
        section (Section): Mutable morphology's section
        linear_density (float): Number of points per micron
    """
    points = section.points
    ids, fractions = _resample_from_linear_density(points, linear_density)
    section.points = _parametric_values(points, ids, fractions)
    section.diameters = _parametric_values(section.diameters, ids, fractions)


def _resample_astrocyte_section(section, linear_density):
    """Resample in-place the section data based on linear density.
    Args:
        section (Section): Mutable morphology's section
        linear_density (float): Number of points per micron
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
        morphio.CellFamily.FAMILY_NEURON: _resample_neuron_section,
        morphio.CellFamily.FAMILY_GLIA: _resample_astrocyte_section
    }[cell_family]


def resample_linear_density(obj, linear_density):
    """Returns a new morphology with new points and point properties,
    the number of which depends on the linear_density. The new
    values are linearly interpolated from the old ones.

    This function is useful for use cases where morphologies have too
    many points and a resampling with fewer points/properties is required.

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
