'''The morph-tool command line launcher'''
import logging
import sys
from pathlib import Path

import click
import dask.bag as dask_bag

import morph_tool
from morph_tool import converter
from morph_tool.utils import iter_morphology_files

logging.basicConfig()
L = logging.getLogger('morph_tool')
L.setLevel(logging.INFO)

REQUIRED_PATH = click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True)


@click.group()
def cli():
    '''The CLI entry point'''


@cli.group()
def convert():
    '''Convert a file format between its different representation.

    A special care has been given to the soma conversion as each file format
    has its own representation of the soma.
    While the soma shape cannot be preserved during the conversion,
    the soma surface should be.

    More information at: https://bbpteam.epfl.ch/project/issues/browse/NSETM-458'''


@cli.command(short_help='Get soma surface as computed by NEURON')
@click.argument('input_file', type=REQUIRED_PATH)
@click.option('--quiet/--no-quiet', default=False)
def soma_surface(input_file, quiet):
    '''Get soma surface as computed by NEURON'''
    # pylint: disable=import-outside-toplevel
    from morph_tool.neuron_surface import get_NEURON_surface
    if quiet:
        L.setLevel(logging.WARNING)

    try:
        click.echo('Soma surface: {}'.format(get_NEURON_surface(input_file)))
    except ImportError as e:
        raise ImportError('''the NEURON module is not installed.
        - if you are on the cluster, you can try: module load nix/hpc/neuron
        - otherwise, get it here: https://github.com/neuronsimulator/nrn and compile it...'''
                          ) from e


@convert.command(short_help='Convert a single morphology')
@click.argument('input_file', type=REQUIRED_PATH)
@click.argument('output_file')
@click.option('--quiet/--no-quiet', default=False)
@click.option('--recenter', is_flag=True,
              help='recenter the morphology based on the center of gravity of the soma')
@click.option('--nrn-order', is_flag=True,
              help='whether to traverse the neuron in the NEURON fashion')
@click.option('--single-point-soma', is_flag=True,
              help='For SWC files only')
def file(input_file, output_file, quiet, recenter, nrn_order, single_point_soma):
    '''Convert a single morphology from/to the following formats: ASC, SWC, H5'''
    if quiet:
        L.setLevel(logging.WARNING)

    converter.convert(input_file, output_file, recenter, nrn_order, single_point_soma)


def _attempt_convert(path, output_dir, extension, recenter, nrn_order, single_point_soma):
    '''Function to be passed to dask.bag.map

    Attempts a conversion and returns the path if it failed
    '''
    try:
        converter.convert(path, Path(output_dir) / (path.stem + '.' + extension),
                          recenter, nrn_order, single_point_soma)
        return None
    except:  # noqa, pylint: disable=bare-except
        return str(path)


@convert.command(short_help='Convert all morphologies in a folder')
@click.argument('input_dir')
@click.argument('output_dir', type=click.Path(exists=True, file_okay=False, writable=True))
@click.option('-ext', '--extension', type=click.Choice(['h5', 'swc', 'asc', 'H5', 'SWC', 'ASC']),
              help='The output file format')
@click.option('--quiet/--no-quiet', default=False)
@click.option('--recenter', is_flag=True,
              help='recenter the morphology based on the center of gravity of the soma')
@click.option('--nrn-order', is_flag=True,
              help='whether to traverse the neuron in the NEURON fashion')
@click.option('--single-point-soma', is_flag=True,
              help='For SWC files only')
@click.option('--ncores', help='The number of cores', default=None)
def folder(input_dir, output_dir, extension, quiet, recenter, nrn_order, single_point_soma, ncores):
    '''Convert all morphologies in the folder and its subfolders'''
    if quiet:
        L.setLevel(logging.WARNING)

    failed_conversions = dask_bag.from_sequence(iter_morphology_files(input_dir),
                                                npartitions=ncores).map(
        _attempt_convert,
        output_dir=output_dir,
        extension=extension,
        recenter=recenter,
        nrn_order=nrn_order,
        single_point_soma=single_point_soma)
    failed_conversions = list(filter(None, failed_conversions))

    if failed_conversions:
        L.warning('The following morphologies could not be converted: %s',
                  failed_conversions)


@cli.command()
@click.argument('morph1')
@click.argument('morph2')
@click.option('--quiet/--no-quiet', default=False)
def diff(morph1, morph2, quiet):
    '''
    Compare two morphology files.

    Return exit code 0 if morphology objects stored in `file1` and `file2`
    are deemed identical by MorphIO; and exit code 1 otherwise.

    Morphologies with different formats can be compared.

    MORPH1 and MORPH2 can be either filenames or morphio.Morpholgy objects
    '''
    if quiet:
        L.setLevel(logging.WARNING)

    result = morph_tool.diff(morph1, morph2)
    if result:
        L.info("Morphologies not identical")
        L.info(result.info)
        sys.exit(1)
