"""Module for morphology simplification."""

import collections
import logging

import numpy as np
import morphio.mut


L = logging.getLogger("morph_tool")


def _squared_distance_points_to_line(points, line_start, line_end):
    """Calculate squared distance of point to line using vector projections.

    The vector projections of the points to the line are calculated and subtracted from the point
    vectors to get the vector rejections, the squared norm of which corresponds to the distance
    squared from the points to the line.

    See: https://en.wikipedia.org/wiki/Vector_projection
    """
    vectors = points - line_start
    line_vector = line_end - line_start

    inv_sq_norm_line = 1.0 / line_vector.dot(line_vector)
    vector_projections = (
        line_vector * vectors.dot(line_vector)[:, np.newaxis] * inv_sq_norm_line
    )

    # get the rejection vectors perpendicular of the projections and get their squared norm
    return np.square(vector_projections - vectors).sum(axis=1)


def _ramer_douglas_peucker(points, epsilon):
    """Simplify polyline determined by 'points' using the Ramer-Douglas-Peucker algorithm.

    Args:
        points (np.ndarray): A sequence of points representing a polyline to be simplified.
        epsilon (float): Distance threshold that determines the max distance between a line segment
            and the original curve.

    Returns:
        A boolean mask corresponding tot the points to keep.
    """
    # keep all points if epsilon is zero
    if np.isclose(epsilon, 0.0):
        return np.ones(len(points), dtype=bool)

    keep_mask = np.zeros(len(points), dtype=bool)

    # first and last points always kept
    keep_mask[[0, -1]] = True

    squared_epsilon = epsilon**2

    beg_index = 0
    end_index = len(points) - 1

    stack = collections.deque([(beg_index, end_index)])

    while stack:

        # start and end index of the current line
        beg_index, end_index = stack.pop()

        # cannot further simplify
        if end_index - beg_index < 3:
            continue

        # calculate the squared distances of all points to the current line
        # the first and last points are not included in the calc because they are always kept
        squared_distances = _squared_distance_points_to_line(
            points=points[beg_index + 1: end_index],
            line_start=points[beg_index],
            line_end=points[end_index],
        )

        # furthest point from current line
        index = np.argmax(squared_distances)

        if squared_distances[index] > squared_epsilon:

            # +1 to account for that first point we didn't include in the distance calculation
            global_index = beg_index + index + 1

            keep_mask[global_index] = True

            # divide the line into two lines and repeat
            stack.append((global_index, end_index))
            stack.append((beg_index, global_index))

    return keep_mask


def simplify_morphology(morph, epsilon):
    """Simplify the sections of a morphology.

    Each section in the morphology is simplified using the Ramer-Douglas-Peucker algorithm.

    The algorithm starts from the first and last points of a curve represented as a list of points
    and recursively divides them into smaller segments.

    At each step it finds the furthest point from the segment. If the point is not within a certain
    threshold distance 'epsilon' from the segment then this point can not be simplified and the
    segment is split into two smaller ones passing through that point, else all the intermediate
    points of this segment are discarded. The algorithm is then applied recursively to each new
    segment.

    The result is a simplified curve that approximates the original with a fewer number of points.
    The distance threshold controls the simplification strength. Higher epsilon values result in
    simpler and less accurate curves with less inflection points (turns). A zero epsilon will result
    in a copy of the initial morphology, whereas a large epsilon will only keep the start and end
    points of each section.

    Args:
        morph: Morphology object, mutable or immutable
        epsilon (float): distance tolerance in microns

    Returns:
        morphio.mut.Morphology: Simplified morphology copy

    Note:
        This algorithm simplifies morphologies, therefore it will change morphometrics such as the
        tortuosity and the path length.
    """
    morph = morphio.mut.Morphology(morph)

    for section in morph.iter():

        points = section.points

        if len(points) == 2:
            continue

        keep_mask = _ramer_douglas_peucker(points, epsilon)

        section.points = points[keep_mask]
        section.diameters = section.diameters[keep_mask]

        if section.perimeters:
            section.perimeters = section.perimeters[keep_mask]

    return morph
