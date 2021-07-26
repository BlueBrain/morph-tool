from pathlib import Path
import itertools as it
import numpy as np

import pytest

import morphio
from morph_tool import diff
from morph_tool.exceptions import MorphToolException
from morph_tool.converter import convert

DATA = Path(__file__).parent / 'data'


def test_convert(tmpdir):
    for in_ext, out_ext in it.product(['asc', 'h5', 'swc'], repeat=2):
        # A simple morphology
        inname = DATA / f'simple.{in_ext}'
        outname = Path(tmpdir, f'test.{out_ext}')
        convert(inname, outname)
        assert not diff(inname, outname)

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
    np.testing.assert_approx_equal(m.soma.surface, 476.0504050847511)


def test_convert_sanitize(tmpdir):
    # needs to have a complex contour soma
    simple = DATA / 'single_child.asc'
    outname = Path(tmpdir, 'single_child.swc')
    with pytest.raises(MorphToolException, match='Use `sanitize` option for converting'):
        convert(simple, outname, single_point_soma=True)

    convert(simple, outname, single_point_soma=True, sanitize=True)
    m = morphio.Morphology(outname)
    assert len(m.sections) == 1
