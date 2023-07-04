"""Test amira_converter module."""
from pathlib import Path
from morph_tool import amira_converter
from morphio.mut import Morphology
from numpy.testing import assert_array_almost_equal


DATA = Path(__file__).parent / "data"


def test_load_amira(tmpdir):
    m_amira = amira_converter.load_amira(DATA / "amira_cell.am")
    m_amira.write(tmpdir / 'amira_cell.asc')
    m_asc = Morphology(DATA / "amira_cell.asc")
    for section_amira, section_asc in zip(m_amira.iter(), m_asc.iter()):
        assert_array_almost_equal(section_amira.points, section_asc.points)
    assert_array_almost_equal(m_amira.soma.points, m_asc.soma.points)
