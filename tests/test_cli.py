import itertools as it
import os
import shutil
from pathlib import Path
from nose.tools import assert_equal, ok_
from click.testing import CliRunner
from morph_tool.cli import cli
from utils import setup_tempdir
from morphio import set_maximum_warnings

DATA = Path(__file__).parent / 'data'

def test_cli():
    runner = CliRunner()
    filename = str(DATA / 'simple.asc')
    result = runner.invoke(cli, ['diff', filename, filename])
    assert_equal(result.exit_code, 0)

    result = runner.invoke(cli, ['diff',
                                 '--quiet',
                                 filename,
                                 str(DATA / 'simple.swc')])
    assert_equal(result.exit_code, 1)

def test_convert():
    set_maximum_warnings(0)
    with setup_tempdir('test-convert-file') as tmp_dir:
        runner = CliRunner()
        filename = Path(DATA, 'simple.asc')
        output_name = Path(tmp_dir, 'simple.h5')
        result = runner.invoke(cli, ['convert', 'file', str(filename), str(output_name)])
        assert_equal(result.exit_code, 0, result.exception)
        ok_(output_name.exists())


    with setup_tempdir('test-convert-folder') as tmp_dir:
        runner = CliRunner()
        result = runner.invoke(cli, ['convert', 'folder', '-ext', 'swc',
                                     str(DATA / 'input-convert'), tmp_dir])
        assert_equal(result.exit_code, 0, result.exc_info)

        n_converted_files = len(list(Path(tmp_dir).rglob('**/*.swc')))

        assert_equal(n_converted_files, 2)
