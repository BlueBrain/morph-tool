'''The morph-tool command line launcher'''
import click


@click.group()
def cli():
    '''The CLI entry point'''
    pass


@cli.command(short_help='Converter')
@click.argument('input_file')
@click.argument('output_file')
def converter(input_file, output_file):
    '''The file converter entry point'''
    from morph_tool.converter import run
    run(input_file, output_file)
