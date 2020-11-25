'''Simplify all neurites using the Ramer-Douglas-Peucker, on a per section basis'''
import logging

import numpy as np

from morphio.mut import Morphology, Section

L = logging.getLogger('morph_tool')

def _dist_line2point(x0, start, end):
    '''distance of x0 from line defined by start, to end
        http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    '''
    diff_start_end = end - start
    return np.divide(np.linalg.norm(np.cross(diff_start_end, start - x0)),
                     np.linalg.norm(diff_start_end))

def _ramer_douglas_peucker(points, epsilon):
    ''''''
    max_dist = 0.0
    index = -1

    for i in range(1, len(points)):
        dist = _dist_line2point(points[i], start=points[0], end=points[-1])
        if max_dist < dist:
            index = i
            max_dist = dist

    if epsilon < max_dist:
        r1 = _ramer_douglas_peucker(points[:index + 1, :], epsilon)
        r2 = _ramer_douglas_peucker(points[index:, :], epsilon)
        return np.vstack((r1[:-1], r2))

    return np.vstack((points[0], points[-1]))


def _points_simplify(points, epsilon):
    '''use Ramer-Douglas-Peucker to simplify the points in a section'''
    simplified = _ramer_douglas_peucker(points, epsilon)

    if np.all(points[0] != simplified[0]):
        L.warning('start points mismatch: %s != %s', points[0], simplified[0])

    if np.all(points[-1] != simplified[-1]):
        L.warning('end points mismatch: %s != %s', points[-1], simplified[-1])

    return simplified


def simplify_neuron(morph, epsilon):
    morph = morph.as_mutable()

    for section in morph.iter():
        section.points = _points_simplify(section.points, epsilon)
        #hack - need to return mask of which points used
        section.diameters = section.diameters[:len(section.points)]
        section.perimeters = section.perimeters[:len(section.points)]

    return morph
