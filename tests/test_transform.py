from pathlib import Path
import pytest

import numpy as np
import numpy.testing as npt

import morphio

import morph_tool.transform as tested


DATA = Path(__file__).parent / 'data'


@pytest.fixture
def morph():
    return morphio.mut.Morphology(DATA / 'simple.asc')


def test_sanity(morph):
    """ Ensure test data is what we expect it to be. """
    npt.assert_almost_equal(
        morph.soma.points,
        [[ 0.1,  0. ,  0. ],
         [ 0.1,  0.1,  0. ],
         [-0.1,  0. ,  0. ],
         [-0.1, -0.1,  0. ]],
    )
    npt.assert_almost_equal(
        morph.sections[0].points[-1],
        [0., 5., 0.]
    )
    npt.assert_almost_equal(
        morph.sections[5].points[-1],
        [-5., -4., 0.]
    )

def test_translate_1(morph):
    """ Translate whole morphology inplace. """
    tested.translate(morph, [1., 2., 3.])
    npt.assert_almost_equal(
        morph.soma.points,
        [[1.1,  2. ,  3. ],
         [1.1,  2.1,  3. ],
         [0.9,  2. ,  3. ],
         [0.9,  1.9,  3. ]],
    )
    npt.assert_almost_equal(
        morph.sections[0].points[-1],
        [1., 7., 3.]
    )
    npt.assert_almost_equal(
        morph.sections[5].points[-1],
        [-4., -2., 3.]
    )

def test_translate_2(morph):
    """ Translate axon inplace. """
    axon = morph.root_sections[1]
    tested.translate(axon, [1., 2., 3.])
    npt.assert_almost_equal(
        axon.children[1].points[-1],
        [-4., -2., 3.]
    )
    # soma should remain intact
    npt.assert_almost_equal(
        morph.soma.points,
        [[ 0.1,  0. ,  0. ],
         [ 0.1,  0.1,  0. ],
         [-0.1,  0. ,  0. ],
         [-0.1, -0.1,  0. ]],
    )
    # other neurites should remain intact
    npt.assert_almost_equal(
        morph.sections[0].points[-1],
        [0., 5., 0.]
    )

def test_translate_3(morph):
    """ Check translation vector. """
    with pytest.raises(ValueError):
        tested.translate(morph, np.identity(4))

def test_rotate_1(morph):
    """ Rotate whole morphology inplace. """
    A = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]  # rotate around Z by 90 degrees
    tested.rotate(morph, A)
    npt.assert_almost_equal(
        morph.soma.points[0],
        [0., -0.1, 0.]
    )
    npt.assert_almost_equal(
        morph.sections[0].points[-1],
        [5., 0., 0.]
    )
    npt.assert_almost_equal(
        morph.sections[5].points[-1],
        [-4., 5., 0.]
    )

def test_rotate_different_origin(morph):
    """ Rotate whole morphology inplace. """
    A = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]  # rotate around Z by 90 degrees

    tested.rotate(morph, A, origin=(1, 1, 1))
    npt.assert_almost_equal(
        morph.soma.points,
        [[ 0. ,  1.9,  0. ],
         [ 0.1,  1.9,  0. ],
         [ 0. ,  2.1,  0. ],
         [-0.1,  2.1,  0. ]],
    )
    npt.assert_almost_equal(
        morph.sections[0].points[-1],
        [5., 2., 0.]
    )
    npt.assert_almost_equal(
        morph.sections[5].points,
        [[-4., 2., 0.],
         [-4., 7., 0.]]
    )

def test_rotate_2(morph):
    """ Check rotation matrix shape. """
    with pytest.raises(ValueError):
        tested.rotate(morph, np.identity(4))

def test_transform_1(morph):
    """ Transform whole morphology inplace. """
    # scale Y coordinate by 3 + shift by (-1, -2, -3)
    A = [[1, 0, 0, -1], [0, 3, 0, -2], [0, 0, 1, -3], [0, 0, 0, 1]]
    tested.transform(morph, A)
    npt.assert_almost_equal(
        morph.soma.points,
        [[-0.9, -2. , -3. ],
         [-0.9, -1.7, -3. ],
         [-1.1, -2. , -3. ],
         [-1.1, -2.3, -3. ]]
    )
    npt.assert_almost_equal(
        morph.sections[0].points[-1],
        [-1., 13., -3.]
    )
    npt.assert_almost_equal(
        morph.sections[5].points[-1],
        [-6., -14., -3.]
    )

def test_transform_2(morph):
    """ Check transform matrix shape. """
    with pytest.raises(ValueError):
        tested.transform(morph, np.identity(3))


def test_align(morph):
    section = morph.section(0)
    tested.align(section, [-1, -1, 0])

    npt.assert_almost_equal(section.points,
                            [[0., 0., 0.], [-3.535534, -3.535534, 0.]])

    child = section.children[0]
    npt.assert_almost_equal(child.points,
                            [[-3.535534, -3.535534, 0.], [0, -7.071068, 0]],
                            decimal=5)

def test_align_already_aligned(morph):
    section = morph.section(1)
    tested.align(section, [-1, 0, 0])

    npt.assert_almost_equal(section.points, [[0., 5., 0.], [-5, 5, 0]])

def test_align_pi_angle(morph):
    section = morph.section(1)
    tested.align(section, [+1, 0, 0])

    npt.assert_almost_equal(section.points, [[0., 5., 0.], [+5, 5, 0]])

def test_align_morphology(morph):
    # Test with apical trunk
    morph = morphio.mut.Morphology(DATA / 'apical_test.h5')
    rotation_mat = tested.align_morphology(morph, [0, 0, 1], method="trunk")
    npt.assert_almost_equal(
        rotation_mat,
        [[ 1,  0,  0],
         [ 0,  0, -1],
         [ 0,  1,  0]]
    )

    # Test with apical root section
    morph = morphio.mut.Morphology(DATA / 'apical_test.h5')
    rotation_mat = tested.align_morphology(morph, [0, 0, 1], method="first_section")
    npt.assert_almost_equal(
        rotation_mat,
        [[ 1,  0,  0],
         [ 0,  0, -1],
         [ 0,  1,  0]]
    )

    # Test with apical root segment
    morph = morphio.mut.Morphology(DATA / 'apical_test.h5')
    rotation_mat = tested.align_morphology(morph, [0, 0, 1], method="first_segment")
    npt.assert_almost_equal(
        rotation_mat,
        [[ 1,  0,  0],
         [ 0,  0, -1],
         [ 0,  1,  0]]
    )

    # Test with whole apical
    morph = morphio.mut.Morphology(DATA / 'apical_test.h5')
    rotation_mat = tested.align_morphology(morph, [0, 0, 1], method="whole")
    npt.assert_almost_equal(
        rotation_mat,
        [[ 9.99998738e-01, -1.12359279e-03, -1.12359354e-03],
         [-1.12359279e-03,  1.26246235e-06, -9.99999404e-01],
         [ 1.12359354e-03,  9.99999404e-01,  2.22044605e-16]]
    )

    # Test with whole apical on non-centered morphology
    morph = morphio.mut.Morphology(DATA / 'apical_test.h5')
    shift = [1.0, 2.0, 3.0]
    tested.translate(morph, shift)
    rotation_mat = tested.align_morphology(morph, [0, 0, 1], method="whole")
    npt.assert_almost_equal(
        rotation_mat,
        [[ 9.99998738e-01, -1.12359279e-03, -1.12359354e-03],
         [-1.12359279e-03,  1.26246235e-06, -9.99999404e-01],
         [ 1.12359354e-03,  9.99999404e-01,  2.22044605e-16]]
    )


    # Test with custom direction
    morph = morphio.mut.Morphology(DATA / 'apical_test.h5')
    rotation_mat = tested.align_morphology(morph, [0, 0, 1], method=[0, 1, 0])
    npt.assert_almost_equal(
        rotation_mat,
        [[ 1,  0,  0],
         [ 0,  0, -1],
         [ 0,  1,  0]]
    )
