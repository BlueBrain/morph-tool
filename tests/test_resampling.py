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
        tested._dispatch_section_function(CellFamily.FAMILY_NEURON),
        tested._resample_neuron_section
    )

    assert isinstance(
        tested._dispatch_section_function(CellFamily.FAMILY_GLIA),
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
