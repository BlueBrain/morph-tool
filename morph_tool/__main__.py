'''The morph-tool command line launcher'''
import logging
import sys

import click

logging.basicConfig()
logger = logging.getLogger('morph_tool')
logger.setLevel(logging.INFO)


@click.group()
def cli():
    '''The CLI entry point'''


@cli.command(short_help='Get soma surface as computed by NEURON')
@click.argument('input_file')
@click.option('--quiet/--no-quiet', default=False)
def soma_surface(input_file, quiet):
    '''Get soma surface as computed by NEURON'''
    from morph_tool.neuron_surface import get_NEURON_surface
    if quiet:
        logger.setLevel(logging.WARNING)

    try:
        click.echo('Soma surface: {}'.format(get_NEURON_surface(input_file)))
    except ImportError:
        raise ImportError('''the NEURON module is not installed.
        - if you are on the cluster, you can try: module load nix/hpc/neuron
        - otherwise, get it here: https://github.com/neuronsimulator/nrn and compile it...''')


@cli.command(short_help='Convert files from/to the following formats: ASC, SWC, H5')
@click.argument('input_file')
@click.argument('output_file')
@click.option('--quiet/--no-quiet', default=False)
def convert(input_file, output_file, quiet):
    '''Convert a file format between its different representation.

    A special care has to be given to the soma conversion as each file format
    as its own representation of the soma.
    While the soma shape cannot be preserved, the soma surface should be.

    More information at: https://bbpteam.epfl.ch/project/issues/browse/NSETM-458'''
    from morph_tool.converter import run

    if quiet:
        logger.setLevel(logging.WARNING)
    run(input_file, output_file)


@cli.command()
@click.argument('file1')
@click.argument('file2')
def diff(file1, file2):
    '''
    Compare two morphology files.

    Return exit code 0 if morphology objects stored in `file1` and `file2`
    are deemed identical by MorphIO; and exit code 1 otherwise.

    Morphologies can be stored in different formats.
    '''
    from morphio import Morphology
    morph1 = Morphology(file1)
    morph2 = Morphology(file2)
    if morph1 != morph2:
        logger.info("Morphologies not identical")
        sys.exit(1)
