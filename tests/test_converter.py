from pathlib import Path
import itertools as it
import numpy as np
import numpy.testing as npt

import pytest

import morphio
from morphio import SomaType
from morphio.mut import Morphology

from morph_tool import converter
from morph_tool import diff
from morph_tool.exceptions import MorphToolException
from morph_tool.converter import convert
from morph_tool.neuron_surface import get_NEURON_surface

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


def test_convert_ensure_NRN_area(tmpdir):
    simple = DATA / 'simple.swc'
    outname = Path(tmpdir, 'test.asc')
    convert(simple, outname, ensure_NRN_area=False)

    npt.assert_almost_equal(get_NEURON_surface(simple), 12.5663, decimal=4)
    npt.assert_almost_equal(get_NEURON_surface(outname), 11.9835, decimal=4)

    convert(simple, outname, ensure_NRN_area=True)
    npt.assert_almost_equal(get_NEURON_surface(outname), 12.59102, decimal=4)


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


def test_convert_swc_cylinder_to_contour(tmpdir):
    # needs to have a complex contour soma
    in_morph = morphio.Morphology("""
      1 1  0  0 0 1. -1
      2 1  -1  0 0 1. 1
      3 1  -1  0 0 1. 2
      4 1  0  0 0 1. 3
      5 3  0  0 0 1.  1
      6 3  0  5 0 1.  5
      7 3 -5  5 0 1.5 6
      8 3  6  5 0 1.5 6
      9 2  0  0 0 1.  1
     10 2  0 -4 0 1.  9
     11 2  6 -4 0 2.  10
     12 2 -5 -4 0 2.  10""", "swc")
    inname = Path(tmpdir, 'test.swc')
    in_morph.as_mutable().write(inname)
    outname = Path(tmpdir, 'test.asc')
    convert(inname, outname)

    m = morphio.Morphology(outname)
    assert 8 == len(m.soma.points)
    assert 8 == len(m.soma.diameters)

    # Test with all equal points
    in_morph = morphio.Morphology("""
      1 1  0  0 0 1. -1
      2 1  0  0 0 1. 1
      3 1  0  0 0 1. 2
      4 1  0  0 0 1. 3
      5 3  0  0 0 1.  1
      6 3  0  5 0 1.  5
      7 3 -5  5 0 1.5 6
      8 3  6  5 0 1.5 6
      9 2  0  0 0 1.  1
     10 2  0 -4 0 1.  9
     11 2  6 -4 0 2.  10
     12 2 -5 -4 0 2.  10""", "swc")
    inname = Path(tmpdir, 'test_all_equal.swc')
    in_morph.as_mutable().write(inname)
    outname = Path(tmpdir, 'test_all_equal.asc')
    with pytest.raises(
        MorphToolException,
        match="Can not convert SOMA_CYLINDERS if all points are equal",
    ):
        convert(inname, outname)


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


def test_all_combinations(tmp_path):
    """Test all conversion combinations."""
    soma_types = [getattr(SomaType, i) for i in dir(SomaType) if not i.startswith("_") and i not in ["name", "value"]]

    # Checking valid types for each extension
    errors = []

    for ext in [".asc", ".h5", ".swc"]:
        for soma_type in soma_types:
            if soma_type == SomaType.SOMA_SINGLE_POINT:
                morph = morphio.Morphology("1 1  0  0 0 1. -1", "swc").as_mutable()
                morph.soma.points = [[0, 0, 0]]
                morph.soma.diameters = [1]
                morph.soma.type = soma_type
            elif soma_type == SomaType.SOMA_CYLINDERS:
                morph = morphio.Morphology("1 1  0  0 0 1. -1", "swc").as_mutable()
                morph.soma.points = [[0, 0, 0], [0, 1, 0]]
                morph.soma.diameters = [1] * len(morph.soma.points)
                morph.soma.type = soma_type
            elif soma_type == SomaType.SOMA_NEUROMORPHO_THREE_POINT_CYLINDERS:
                morph = morphio.Morphology("1 1  0  0 0 1. -1", "swc").as_mutable()
                morph.soma.points = [[0, 0, 0], [0, 1, 0], [0, -1, 0]]
                morph.soma.diameters = [2] * len(morph.soma.points)
                morph.soma.type = soma_type
            elif soma_type == SomaType.SOMA_SIMPLE_CONTOUR:
                morph = morphio.Morphology("((Dendrite))", "asc").as_mutable()
                morph.soma.points = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0]]
                morph.soma.diameters = [0] * len(morph.soma.points)
                morph.soma.type = soma_type
            elif soma_type == SomaType.SOMA_UNDEFINED:
                morph = morphio.Morphology("((Dendrite))", "asc").as_mutable()
                morph.soma.points = [[0, 0, 0], [0, 1, 0]]
                morph.soma.diameters = [1] * len(morph.soma.points)
                morph.soma.type = soma_type
            else:
                raise ValueError(f"Unknown soma type: {soma_type}")

            output_file = tmp_path / f"morph_{soma_type.name}{ext}"

            # Convert the morphology
            convert(morph, output_file)

            # Check that the morphology can be read
            Morphology(output_file)
