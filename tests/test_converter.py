import os
from itertools import product

from nose.tools import ok_

from morph_tool import diff
from morph_tool.converter import run
from utils import setup_tempdir

PATH = os.path.join(os.path.dirname(__file__), 'data')


def test_convert():
    with setup_tempdir('test-convert') as tmp_dir:
        for in_ext, out_ext in product(['asc', 'h5', 'swc'], repeat=2):
            # A simple morphology
            inname = os.path.join(PATH, 'simple.' + in_ext)
            outname = os.path.join(tmp_dir, 'test.' + out_ext)
            run(inname, outname)
            ok_(not diff(inname, outname))

            # A more complex one
            inname = os.path.join(PATH, 'tkb061126a4_ch0_cc2_h_zk_60x_1.' + in_ext)
            outname = os.path.join(tmp_dir, 'test.' + out_ext)
            run(inname, outname)
            diff_result = diff(inname, outname, rtol=1e-5, atol=1e-5)
            ok_(not bool(diff_result), diff_result.info)
