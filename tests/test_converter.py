from pathlib import Path
import itertools as it
import numpy as np
import numpy.testing as npt

import pytest

import morphio
from morph_tool import converter
from morph_tool import diff
from morph_tool.exceptions import MorphToolException
from morph_tool.converter import convert

DATA = Path(__file__).parent / 'data'


def test_convert(tmpdir):
    for in_ext, out_ext in (('asc', 'h5'),
                            ('asc', 'swc'),
                            ('h5', 'asc'),
                            ('h5', 'swc')):
        # A simple morphology
        inname = DATA / f'simple.{in_ext}'
        outname = Path(tmpdir, f'test.{out_ext}')
        convert(inname, outname, single_point_soma=(out_ext == 'swc'))
        assert not diff(inname, outname, skip_perimeters=in_ext == 'h5')

        # A more complex one
        inname = DATA / f'tkb061126a4_ch0_cc2_h_zk_60x_1.{in_ext}'
        outname = Path(tmpdir, f'test.{out_ext}')
        convert(inname, outname)
        diff_result = diff(inname, outname, rtol=1e-5, atol=1e-5)
        assert not bool(diff_result), diff_result.info


def test_convert_recenter(tmpdir):
    simple = DATA / 'simple.swc'
    outname = Path(tmpdir, 'test.swc')
    convert(simple, outname, recenter=True)
    assert not diff(simple, outname)  #simple.swc is already centered

    mut = morphio.Morphology(simple).as_mutable()
    mut.soma.points = [[1, 1, 1], ]
    inname = Path(tmpdir, 'moved.swc')
    mut.write(inname)

    convert(inname, outname, recenter=True)
    simple = morphio.Morphology(simple)
    centered_morph = morphio.Morphology(outname)
    assert np.all((simple.points - centered_morph.points) == 1)


def test_convert_swc_contour_to_sphere(tmpdir):
    # needs to have a complex contour soma
    simple = DATA / 'tkb061126a4_ch0_cc2_h_zk_60x_1.asc'
    outname = Path(tmpdir, 'test.swc')
    convert(simple, outname, single_point_soma=True)

    m = morphio.Morphology(outname)
    assert 1 == len(m.soma.points)
    assert 1 == len(m.soma.diameters)

    #value dumped from NEURON: h.area(0.5, sec=icell.soma[0])
    npt.assert_approx_equal(m.soma.surface, 476.0504050847511)


def test_convert_sanitize(tmpdir):
    # needs to have a complex contour soma
    simple = DATA / 'single_child.asc'
    outname = Path(tmpdir, 'single_child.swc')
    with pytest.raises(MorphToolException, match='Use `sanitize` option for converting'):
        convert(simple, outname, single_point_soma=True)

    convert(simple, outname, single_point_soma=True, sanitize=True)
    m = morphio.Morphology(outname)
    assert len(m.sections) == 1


def test__create_contour(tmpdir):
    for point_count in range(3, 20):
        points, diameters = converter._create_contour(radius=10, point_count=point_count, line_width=0.1)
        assert len(points) == point_count
        assert len(diameters) == point_count
        npt.assert_allclose(points[:, 2], [0]*point_count)
        npt.assert_allclose(diameters, [0.1]*point_count)
        assert len(np.unique(np.around(points, decimals=4), axis=0)) == point_count
