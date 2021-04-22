'''A utility module to get the surface as computed by neuron'''


def get_NEURON_surface(path):
    '''Return the soma surface computed by the NEURON simulator. Maybe.'''
    try:
        # pylint: disable=import-outside-toplevel
        from bluepyopt.ephys import simulators, morphologies, models  # pylint: disable=import-error
    except ImportError as e:
        raise ImportError(
            'bluepyopt not installed; please use `pip install morph-tool[nrn]`') from e

    SIM = simulators.NrnSimulator()
    HOC = SIM.neuron.h

    icell = models.CellModel('model', path, [], [])
    morph = morphologies.NrnFileMorphology(
        path, do_replace_axon=False, do_set_nseg=True)
    morph.instantiate(SIM, icell)
    surface = sum(HOC.area(0.5, sec=soma_sec) for soma_sec in icell.soma)  # noqa pylint: disable=no-member
    return surface
