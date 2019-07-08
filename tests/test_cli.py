import shutil
from os.path import dirname, join as joinp
from nose.tools import assert_equal
from click.testing import CliRunner
from morph_tool.cli import cli


PATH = joinp(dirname(__file__), 'data')


def test_cli():
    runner = CliRunner()
    filename = joinp(PATH, 'simple.asc')
    result = runner.invoke(cli, ['diff', filename, filename])
    assert_equal(result.exit_code, 0)

    result = runner.invoke(cli, ['diff', '--quiet', filename, joinp(PATH, 'simple.swc')])
    assert_equal(result.exit_code, 1)
