import itertools as it
import os
import shutil
from nose.tools import assert_equal
from click.testing import CliRunner
from morph_tool.cli import cli
from utils import setup_tempdir

PATH = os.path.join(os.path.dirname(__file__), 'data')


def test_cli():
    runner = CliRunner()
    filename = os.path.join(PATH, 'simple.asc')
    result = runner.invoke(cli, ['diff', filename, filename])
    assert_equal(result.exit_code, 0)

    result = runner.invoke(cli, ['diff',
                                 '--quiet',
                                 filename,
                                 os.path.join(PATH, 'simple.swc')])
    assert_equal(result.exit_code, 1)
