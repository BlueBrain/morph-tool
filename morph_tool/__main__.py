'''The morph-tool command line launcher'''
import click


@click.group()
def cli():
    '''The CLI entry point'''
    pass


@cli.command(short_help='Convert files from/to the following formats: ASC, SWC, H5')
@click.argument('input_file')
@click.argument('output_file')
def convert(input_file, output_file):
    '''The file converter entry point'''
    from morph_tool.converter import run
    run(input_file, output_file)
