import os

import nose.tools as nt

import numpy as np
import numpy.testing as npt

import morphio

import morph_tool.transform as tested


TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def _test_data_path(filename):
    return os.path.join(TEST_DIR, 'data', filename)


class TestTransform:
    def setup(self):
        self.morph = morphio.mut.Morphology(_test_data_path('simple.asc'))

    def test_sanity(self):
        """ Ensure test data is what we expect it to be. """
        npt.assert_almost_equal(
            self.morph.soma.points[0],
            [0., 0., 0.]
        )
        npt.assert_almost_equal(
            self.morph.sections[0].points[-1],
            [0., 5., 0.]
        )
        npt.assert_almost_equal(
            self.morph.sections[5].points[-1],
            [-5., -4., 0.]
        )

    def test_translate_1(self):
        """ Translate whole morphology inplace. """
        tested.translate(self.morph, [1., 2., 3.])
        npt.assert_almost_equal(
            self.morph.soma.points[0],
            [1., 2., 3.]
        )
        npt.assert_almost_equal(
            self.morph.sections[0].points[-1],
            [1., 7., 3.]
        )
        npt.assert_almost_equal(
            self.morph.sections[5].points[-1],
            [-4., -2., 3.]
        )

    def test_translate_2(self):
        """ Translate axon inplace. """
        axon = self.morph.root_sections[1]
        tested.translate(axon, [1., 2., 3.])
        npt.assert_almost_equal(
            axon.children[1].points[-1],
            [-4., -2., 3.]
        )
        # soma should remain intact
        npt.assert_almost_equal(
            self.morph.soma.points[0],
            [0., 0., 0.]
        )
        # other neurites should remain intact
        npt.assert_almost_equal(
            self.morph.sections[0].points[-1],
            [0., 5., 0.]
        )

    def test_translate_3(self):
        """ Check translation vector. """
        nt.assert_raises(
            ValueError,
            tested.translate, self.morph, np.identity(4)
        )

    def test_rotate_1(self):
        """ Rotate whole morphology inplace. """
        A = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]  # rotate around Z by 90 degrees
        tested.rotate(self.morph, A)
        npt.assert_almost_equal(
            self.morph.soma.points[0],
            [0., 0., 0.]
        )
        npt.assert_almost_equal(
            self.morph.sections[0].points[-1],
            [5., 0., 0.]
        )
        npt.assert_almost_equal(
            self.morph.sections[5].points[-1],
            [-4., 5., 0.]
        )

    def test_rotate_different_origin(self):
        """ Rotate whole morphology inplace. """
        A = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]  # rotate around Z by 90 degrees

        tested.rotate(self.morph, A, origin=(1, 1, 1))
        npt.assert_almost_equal(
            self.morph.soma.points[0],
            [0., 2., 0.]
        )
        npt.assert_almost_equal(
            self.morph.sections[0].points[-1],
            [5., 2., 0.]
        )
        npt.assert_almost_equal(
            self.morph.sections[5].points,
            [[-4., 2., 0.],
             [-4., 7., 0.]]
        )

    def test_rotate_2(self):
        """ Check rotation matrix shape. """
        nt.assert_raises(
            ValueError,
            tested.rotate, self.morph, np.identity(4)
        )

    def test_transform_1(self):
        """ Transform whole morphology inplace. """
        # scale Y coordinate by 3 + shift by (-1, -2, -3)
        A = [[1, 0, 0, -1], [0, 3, 0, -2], [0, 0, 1, -3], [0, 0, 0, 1]]
        tested.transform(self.morph, A)
        npt.assert_almost_equal(
            self.morph.soma.points[0],
            [-1., -2., -3.]
        )
        npt.assert_almost_equal(
            self.morph.sections[0].points[-1],
            [-1., 13., -3.]
        )
        npt.assert_almost_equal(
            self.morph.sections[5].points[-1],
            [-6., -14., -3.]
        )

    def test_transform_2(self):
        """ Check transform matrix shape. """
        nt.assert_raises(
            ValueError,
            tested.transform, self.morph, np.identity(3)
        )


    def test_align(self):
        section = self.morph.section(0)
        tested.align(section, [-1, -1, 0])

        npt.assert_almost_equal(section.points,
                                [[0., 0., 0.], [-3.535534, -3.535534, 0.]])

        child = section.children[0]
        npt.assert_almost_equal(child.points,
                                [[-3.535534, -3.535534, 0.], [0, -7.071068, 0]],
                                decimal=5)

    def test_align_already_aligned(self):
        section = self.morph.section(1)
        tested.align(section, [-1, 0, 0])

        npt.assert_almost_equal(section.points, [[0., 5., 0.], [-5, 5, 0]])

    def test_align_pi_angle(self):
        section = self.morph.section(1)
        tested.align(section, [+1, 0, 0])

        npt.assert_almost_equal(section.points, [[0., 5., 0.], [+5, 5, 0]])
