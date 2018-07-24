'''The morph-tool command line launcher'''
import logging

import click

logging.basicConfig()
logger = logging.getLogger('morph_tool')
logger.setLevel(logging.INFO)


@click.group()
def cli():
    '''The CLI entry point'''
    pass


@cli.command(short_help='Get soma surface as computed by NEURON')
@click.argument('input_file')
def soma_surface(input_file):
    '''Get soma surface as computed by NEURON'''
    from morph_tool.neuron_surface import get_NEURON_surface
    try:
        click.echo('Soma surface: {}'.format(get_NEURON_surface(input_file)))
    except ImportError:
        raise ImportError('''the NEURON module is not installed.
        - if you are on the cluster, you can try: module load nix/hpc/neuron
        - otherwise, get it here: https://github.com/neuronsimulator/nrn and compile it...''')


@cli.command(short_help='Convert files from/to the following formats: ASC, SWC, H5')
@click.argument('input_file')
@click.argument('output_file')
def convert(input_file, output_file):
    '''Convert a file format between its different representation.

    A special care has to be given to the soma conversion as each file format
    as its own representation of the soma.
    While the soma shape cannot be preserved, the soma surface should be.

    More information at: https://bbpteam.epfl.ch/project/issues/browse/NSETM-458'''
    from morph_tool.converter import run

    run(input_file, output_file)
