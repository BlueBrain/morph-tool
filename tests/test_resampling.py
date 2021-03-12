from copy import deepcopy
from pathlib import Path
from mock import Mock
import numpy as np
from numpy import testing as npt
from morphio import Morphology
from morph_tool import resampling as tested
from morph_tool.morphio_diff import diff


DATA_DIR = Path(__file__).resolve().parent / 'data'


def test_accumulated_path_lengths():

    points = np.array([[0., 0., 0.], [1., 1., 1.], [2., 2., 2.]])
    segment_length = np.linalg.norm(points[1] - points[0])
    total_length = np.linalg.norm(points[1:] - points[:-1], axis=1).sum()

    path_lengths = tested._accumulated_path_lengths(points)

    npt.assert_allclose(path_lengths, [0.0, segment_length, 2.0 * segment_length])

    # last path length should be equal to the total length
    npt.assert_allclose(path_lengths[-1], total_length)


def test_resample_from_linear_density():
    """
                 | ---------- | ---------- |
    ids   :      0            1            2
    points: [0., 0., 0.] [1., 0., 0.] [2., 0., 0.]
    """
    points = np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]])

    ids, fractions = tested._resample_from_linear_density(points, 2.0)

    # If we don't take into account the start and end we expect a point
    # at every 0.5 um. Therefore:
    # [0.5, 0., 0.], [1., 0., 0.], [1.5, 0., 0.]

    expected_ids = [0, 0, 1]
    expected_fractions = [0.5, 1.0, 0.5]

    # so that p' = points[ids] + fractions * (points[ids + 1] - points[ids])
    # p0 = [0., 0., 0.] + 0.5 * ([1., 0., 0.] - [0., 0., 0.]) = [0.5, 0., 0.]
    # p1 = [0., 0., 0.] + 1.0 * ([1., 0., 0.] - [0., 0., 0.]) = [1., 0., 0.]
    # p2 = [1., 0., 0.] + 0.5 * ([2., 0., 0.] - [1., 0., 0.]) = [1.5, 0., 0.]

    npt.assert_array_equal(ids, expected_ids)
    npt.assert_allclose(fractions, expected_fractions)


def test_resample_parametric_values():

    ids = np.array([0, 0, 1])
    fractions = np.array([0.5, 1.0, 0.5])

    values = np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]])

    # Using and array of 3D points
    new_values = tested._parametric_values(values, ids, fractions)

    npt.assert_allclose(new_values,
        [[0.0, 0., 0.],
         [0.5, 0., 0.],
         [1.0, 0., 0.],
         [1.5, 0., 0.],
         [2.0, 0., 0.]]
    )

    # using an array of 1D values
    values = np.array([1., 2., 3.])
    new_values = tested._parametric_values(values, ids, fractions)

    npt.assert_allclose(new_values, [1.0, 1.5, 2.0, 2.5, 3.0])


def test_resample_neuron_section_strip():
    """Check that only the first and last points and diameters remain
    if linear_density is zero.
    """
    section = Mock(
        points=np.array([[0., 0., 0.], [1., 1., 1.], [2., 2., 2.]]),
        diameters=np.array([3., 2., 1.])
    )

    tested._resample_neuron_section(section, linear_density=0.)

    npt.assert_allclose(section.points, [[0., 0., 0.], [2., 2., 2.]])
    npt.assert_allclose(section.diameters, [3., 1.])


def test_resample_astrocyte_section_strip():
    """Check that only the first and last points and diameters remain
    if linear_density is zero.
    """
    section = Mock(
        points=np.array([[0., 0., 0.], [1., 1., 1.], [2., 2., 2.]]),
        diameters=np.array([3., 2., 1.]),
        perimeters=np.array([1., 2., 3.])
    )

    tested._resample_astrocyte_section(section, linear_density=0.)

    npt.assert_allclose(section.points, [[0., 0., 0.], [2., 2., 2.]])
    npt.assert_allclose(section.diameters, [3., 1.])
    npt.assert_allclose(section.perimeters, [1., 3.])


def test_resample_neuron_section__identity():
    """If we have points 1um apart with a density of 1 points per um
    then the result should the same as the input.
    """
    section = Mock(
        points=np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
        diameters=np.array([3., 2., 1.])
    )

    tested._resample_neuron_section(section, linear_density=1.)

    npt.assert_allclose(
        section.points, [[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]])
    npt.assert_allclose(section.diameters, [3., 2., 1.])


def test_resample_astrocyte_section__identity():
    """If we have points 1um apart with a density of 1 points per um
    then the result should the same as the input.
    """
    section = Mock(
        points=np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
        diameters=np.array([3., 2., 1.]),
        perimeters=np.array([1., 2., 3.])
    )

    tested._resample_astrocyte_section(section, linear_density=1.)

    npt.assert_allclose(
        section.points, [[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]])
    npt.assert_allclose(section.diameters, [3., 2., 1.])
    npt.assert_allclose(section.perimeters, [1., 2., 3.])


def test_resample_neuron_section__double_density():
    """If we have points 1um apart with a density of 1 points per um
    then the result should the same as the input.
    """

    section = Mock(
        points=np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
        diameters=np.array([3., 2., 1.])
    )

    tested._resample_neuron_section(section, linear_density=2.)

    npt.assert_allclose(
        section.points,
        [[0.0, 0., 0.],
         [0.5, 0., 0.],
         [1.0, 0., 0.],
         [1.5, 0., 0.],
         [2.0, 0., 0.]])
    npt.assert_allclose(section.diameters, [3., 2.5, 2., 1.5, 1.])

    # let's do a full cycle now!
    tested._resample_neuron_section(section, linear_density=1.)

    npt.assert_allclose(
        section.points,
        [[0.0, 0., 0.],
         [1.0, 0., 0.],
         [2.0, 0., 0.]])
    npt.assert_allclose(section.diameters, [3., 2., 1.])


def test_resample_astrocyte_section__double_density():
    """If we have points 1um apart with a density of 1 points per um
    then the result should the same as the input.
    """

    section = Mock(
        points=np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
        diameters=np.array([3., 2., 1.]),
        perimeters=np.array([1., 2., 3.])
    )

    tested._resample_astrocyte_section(section, linear_density=2.)

    npt.assert_allclose(
        section.points,
        [[0.0, 0., 0.],
         [0.5, 0., 0.],
         [1.0, 0., 0.],
         [1.5, 0., 0.],
         [2.0, 0., 0.]])
    npt.assert_allclose(section.diameters, [3., 2.5, 2., 1.5, 1.])
    npt.assert_allclose(section.perimeters, [1., 1.5, 2., 2.5, 3.])

    # let's do a full cycle now!
    tested._resample_astrocyte_section(section, linear_density=1.)

    npt.assert_allclose(
        section.points,
        [[0.0, 0., 0.],
         [1.0, 0., 0.],
         [2.0, 0., 0.]])
    npt.assert_allclose(section.diameters, [3., 2., 1.])
    npt.assert_allclose(section.perimeters, [1., 2., 3.])


def _test_dispatch_section_function():

    from morphio import CellFamily

    assert isinstance(
        tested._dispatch_section_function(CellFamily.NEURON),
        tested._resample_neuron_section
    )

    assert isinstance(
        tested._dispatch_section_function(CellFamily.GLIA),
        tested._resample_astrocyte_section
    )

def _assert_allclose_first_last(arr1, arr2):
    npt.assert_allclose(arr1[[0, -1]], arr2[[0, -1]])


def test_resample_linear_density__astrocyte():
    """Run the function and make sure that it doesn't mutate the input cell
    """
    obj = Morphology(DATA_DIR / 'neuron.h5')
    new_obj = tested.resample_linear_density(obj, linear_density=0.).as_immutable()

    # densty 0 result to only first and last data points
    assert len(obj.points) > len(new_obj.points)
    assert len(obj.diameters) > len(new_obj.diameters)
    assert len(obj.perimeters) > len(new_obj.perimeters)

    for s1, s2 in zip(obj.iter(), new_obj.iter()):

        s2_points = s2.points
        s2_diameters = s2.diameters

        assert len(s2_points) == len(s2_diameters) == 2

        _assert_allclose_first_last(s1.points, s2_points)
        _assert_allclose_first_last(s1.diameters, s2_diameters)


def test_resample_linear_density__astrocyte():
    """Run the function and make sure that it doesn't mutate the input cell
    """
    obj = Morphology(DATA_DIR / 'astrocyte.h5')
    new_obj = tested.resample_linear_density(obj, linear_density=0.).as_immutable()

    # densty 0 result to only first and last data points
    assert len(obj.points) > len(new_obj.points)
    assert len(obj.diameters) > len(new_obj.diameters)
    assert len(obj.perimeters) > len(new_obj.perimeters)

    for s1, s2 in zip(obj.iter(), new_obj.iter()):

        s2_points = s2.points
        s2_diameters = s2.diameters
        s2_perimeters = s2.perimeters

        assert len(s2_points) == len(s2_diameters) == len(s2_perimeters) == 2

        _assert_allclose_first_last(s1.points, s2_points)
        _assert_allclose_first_last(s1.diameters, s2_diameters)
        _assert_allclose_first_last(s1.perimeters, s2_perimeters)


def test_resample_from_linear_density__numerical_innacurracy():
    """
    Previous implementation was using arange to calculate the new paths,
    which led to the accumulation of a numerical error. As a result the
    last path was numerically smaller that the total length in the resampling
    function which led in the generation of an extra point.
    """
    points = np.array([
        [-1.474797606, -12.873973846, 1.465831757],
        [-1.268972993, -13.383681297, 1.386174798],
        [-1.207281470, -13.681849480, 1.339448690],
        [-0.772920549, -14.285424232, 1.751677036],
        [-0.640972078, -14.575642586, 1.706354022],
        [-0.157932609, -15.903386116, 1.498738170],
        [0.108525813, -16.787338257, 1.856412053],
        [0.603999615, -17.670700073, 1.718663454],
        [1.443708062, -19.519208908, 1.925997615],
        [1.577559590, -19.740325928, 1.891538978],
        [2.688630342, -22.260030746, 2.004102230],
        [3.322759628, -23.588979721, 2.292642832],
        [4.462557793, -25.080135345, 2.556796074],
        [4.665114880, -25.744903564, 2.948777676]
    ])

    total = np.linalg.norm(points[:1] - points[:-1], axis=1).sum()

    ids, fractions = tested._resample_from_linear_density(points, linear_density=1.)

    new_points = tested._parametric_values(points, ids, fractions)

    assert not np.allclose(new_points[-1], new_points[-2])
