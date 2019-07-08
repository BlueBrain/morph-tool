import shutil
from os.path import dirname, join as joinp
from nose.tools import assert_equal
from click.testing import CliRunner
from morph_tool import diff
from morphio import LogLevel


PATH = joinp(dirname(__file__), 'data')


def test_diff():
    filename = joinp(PATH, 'simple.asc')

    assert_equal(diff(filename, filename), False)

    assert_equal(diff(filename, joinp(PATH, 'simple.swc'),
                      verbose_level=LogLevel.error), True)
