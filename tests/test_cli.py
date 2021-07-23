from pathlib import Path
from click.testing import CliRunner
from morph_tool.cli import cli
from morphio import set_maximum_warnings

DATA = Path(__file__).parent / 'data'

def test_cli():
    runner = CliRunner()
    filename = str(DATA / 'simple.asc')
    result = runner.invoke(cli, ['diff', filename, filename])
    assert result.exit_code == 0

    result = runner.invoke(cli, ['diff',
                                 '--quiet',
                                 filename,
                                 str(DATA / 'simple.swc')])
    assert result.exit_code == 1


def test_convert_file(tmpdir):
    set_maximum_warnings(0)
    runner = CliRunner()
    filename = Path(DATA, 'simple.asc')
    output_name = Path(tmpdir, 'simple.h5')
    result = runner.invoke(cli, ['convert', 'file', str(filename), str(output_name)])
    assert result.exit_code == 0, result.exception
    assert output_name.exists()


def test_convert_folder(tmpdir):
    runner = CliRunner()
    result = runner.invoke(cli, ['convert', 'folder', '-ext', 'swc',
                                 str(DATA / 'input-convert'), str(tmpdir)])
    assert result.exit_code == 0, result.exc_info

    n_converted_files = len(list(Path(tmpdir).rglob('**/*.swc')))

    assert n_converted_files == 2
