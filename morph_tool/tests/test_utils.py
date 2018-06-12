import math
import sys

from bluepyopt.ephys import simulators, morphologies, models


def get_NEURON_surface(path):
    '''Return the soma surface computed by the NEURON simulator'''
    sim = simulators.NrnSimulator()
    h = sim.neuron.h
    icell = models.CellModel('model', path, [], [])
    morph = morphologies.NrnFileMorphology(path, do_replace_axon=False, do_set_nseg=True)
    morph.instantiate(sim, icell)

    assert len(icell.soma) == 1

    soma_sec = icell.soma[0]
    return h.area(0.5, sec=soma_sec)
