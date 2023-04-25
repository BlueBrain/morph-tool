from pathlib import Path

import numpy as np
import pytest
from morph_tool import simplify as test_module
from morphio import Morphology
from numpy import testing as npt

DATA_DIR = Path(__file__).resolve().parent / "data"


@pytest.mark.parametrize(
    "point, lbeg, lend, expected_distance",
    [
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0),  # same point
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.0),  # same point
        ([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.0),  # on line
        ([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [-1.0, -1.0, -1.0], 0.0),  # on line
        ([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [-1.0, -1.0, -1.0], 0.0),  # on line
        ([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),  # perpendicular
        ([0.0, -1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),  # perpendicular
        ([4.0, 6.0, 0.0], [-7.0, 9.0, 0.0], [10.0, 9.0, 0.0], 9.0),  # perpendicular
    ],
)
def test_squared_distance_points_to_vectors(point, lbeg, lend, expected_distance):

    res = test_module._squared_distance_points_to_line(
        points=np.array([point]),
        line_start=np.array(lbeg),
        line_end=np.array(lend),
    )
    npt.assert_almost_equal(res, expected_distance)


def test_ramer_douglas_peucker__zero_epsilon():

    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0],
        ]
    )

    result = points[test_module._ramer_douglas_peucker(points.copy(), epsilon=0.0)]

    npt.assert_allclose(result, points)


def test_ramer_douglas_peucker__large_epsilon():

    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0],
        ]
    )

    result = points[test_module._ramer_douglas_peucker(points.copy(), epsilon=1000.0)]

    npt.assert_allclose(result, points[[0, -1]])


def test_rame_douglas_peucker__sine_function():

    dt = np.pi / 4.0
    t = np.arange(0.0, 2.0 * np.pi + dt, dt)

    points = np.zeros((9, 3), dtype=float)

    points[:, 0] = t
    points[:, 1] = np.sin(t)

    result = points[test_module._ramer_douglas_peucker(points.copy(), epsilon=0.5)]

    expected_points = points[[0, 2, 6, 8]]

    npt.assert_allclose(result, expected_points)


def test_rame_douglas_peucker__expcos_function():

    dt = np.pi / 4.0

    t = np.arange(0.0, 4.0 * np.pi + dt, dt)
    f = lambda t: np.exp(-0.1 * t) * np.sin(t)

    points = np.zeros((len(t), 3), dtype=float)
    points[:, 0] = t
    points[:, 1] = f(t)

    result = points[test_module._ramer_douglas_peucker(points.copy(), epsilon=0.3)]

    expected_points = points[[0, 2, 6, 10, 14, 16]]

    npt.assert_allclose(result, expected_points)


def test_simplify_morphology():

    obj = Morphology(DATA_DIR / "neuron.asc")

    result = test_module.simplify_morphology(obj, 10.0).as_immutable()

    assert len(obj.points) > len(result.points)
