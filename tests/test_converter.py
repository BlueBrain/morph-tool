import numpy as np

import itertools as it
import os

from nose.tools import ok_, assert_equal

import morphio
from morph_tool import diff
from morph_tool.converter import convert
from utils import setup_tempdir

DATA = os.path.join(os.path.dirname(__file__), 'data')


def test_convert():
    with setup_tempdir('test-convert') as tmp_dir:
        for in_ext, out_ext in it.product(['asc', 'h5', 'swc'], repeat=2):
            # A simple morphology
            inname = os.path.join(DATA, 'simple.' + in_ext)
            outname = os.path.join(tmp_dir, 'test.' + out_ext)
            convert(inname, outname)
            ok_(not diff(inname, outname))

            # A more complex one
            inname = os.path.join(DATA, 'tkb061126a4_ch0_cc2_h_zk_60x_1.' + in_ext)
            outname = os.path.join(tmp_dir, 'test.' + out_ext)
            convert(inname, outname)
            diff_result = diff(inname, outname, rtol=1e-5, atol=1e-5)
            ok_(not bool(diff_result), diff_result.info)


def test_convert_recenter():
    with setup_tempdir('test-convert_recenter') as tmp_dir:
        simple = os.path.join(DATA, 'simple.swc')
        outname = os.path.join(tmp_dir, 'test.swc')
        convert(simple, outname, recenter=True)
        ok_(not diff(simple, outname))  #simple.swc is already centered

        mut = morphio.Morphology(simple).as_mutable()
        mut.soma.points = [[1, 1, 1], ]
        inname = os.path.join(tmp_dir, 'moved.swc')
        mut.write(inname)

        convert(inname, outname, recenter=True)
        simple = morphio.Morphology(simple)
        centered_morph = morphio.Morphology(outname)
        ok_(np.all((simple.points - centered_morph.points) == 1))

def test_convert_swc_contour_to_sphere():
    with setup_tempdir('test_convert_swc_contour_to_sphere') as tmp_dir:
        # needs to have a complex contour soma
        simple = os.path.join(DATA, 'tkb061126a4_ch0_cc2_h_zk_60x_1.asc')
        outname = os.path.join(tmp_dir, 'test.swc')
        convert(simple, outname, single_point_soma=True)

        m = morphio.Morphology(outname)
        assert_equal(1, len(m.soma.points))
        assert_equal(1, len(m.soma.diameters))

        #value dumped from NEURON: h.area(0.5, sec=icell.soma[0])
        np.testing.assert_approx_equal(m.soma.surface, 476.0504050847511)
