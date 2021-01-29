The morphdb module introduces two classes:

MorphInfo
=========

A class containing all the information that can be found in the neurondb
regarding a single morphology.

.. code:: python


    info = MorphInfo(name='name', mtype='mtype', 'layer'='layer')
    info.use_dendrite = False
    info.axon_inputs = ['morph1', 'morph2']

On top of the three canonical arguments `name`, `mtype` and `layer`,
`MorphInfo` accepts many keyword arguments:

.. code:: python

   info = MorphInfo(name='name', mtype='mtype', 'layer'='layer',
                    use_axon=False,
                    use_dendrites=False,
                    axon_repair=False,
                    dendrite_repair=False,
                    basal_dendrite_repair=False,
                    tuft_dendrite_repair=False,
                    oblique_dendrite_repair=False,
                    unravel=False,
                    use_for_stats=False,
                    axon_inputs=['morph_1', 'morph_2'],
                    path=Path('/path/to/the/morph'),
                    label='any-arbitrary-label')


MorphDB
=======

A class containing information regarding collections of morphologies.
The class has a single attribute ``self.df`` to expose this information
as a Pandas dataframe.

Constructors:
~~~~~~~~~~~~~

A `MorphDB` object can be built using a sequence of `MorphInfo` objects:

.. code:: python

    db1 = MorphDB([MorphInfo(name=morph1, mtype='L1_DAC', layer=1),
                   MorphInfo(name=morph2, mtype='L2_DAC', layer=1),
                   MorphInfo(name=morph2, mtype='L3_DAC', layer=1)])


.. note::

   Multiple MorphInfo with the same `name` are supported


Or from `neuronDB.xml` or `neuronDB.dat` file:

.. code:: python

          db2 = MorphDB.from_neurondb('neuronDB.xml')
          db2 = MorphDB.from_neurondb('neuronDB.dat')

          # By default, from_neurondb look for morphologies in the folder
          # containing the neuronDB file but this can be overwritten
          db2 = MorphDB.from_neurondb('neuronDB.dat', morphology_folder='/path/to/another/folder')

Or from a standard folder by also passing a list of mtypes

.. code:: python

    db3 = MorphDB.from_folder(my_folder,
                              [(morph1, 'L1_DAC'),
                               (morph2, 'L2_TPC'),
                               (morph2, 'L3_TPC')])

.. note::

   Multiple MorphInfo with the same `name` are supported

Adding new data
~~~~~~~~~~~~~~~

Data can be added through the the ``+`` and ``+=`` operators.

Example:

.. code:: python

    all_morphs = MorphDB.from_neurondb('/gpfs/.../unravelled/neuronDB.xml', label='unravelled')
    all_morphs += MorphDB.from_neurondb('/gpfs/.../repaired/neuronDB.xml', label='repaired')

    more_morphs = MorphDB([MorphInfo(name='name', mtype='mtype', 'layer'='layer', label='extra-morphs'),
                           MorphInfo(name='name', mtype='mtype', 'layer'='layer', label='extra-morphs'),
                           MorphInfo(name='name', mtype='mtype', 'layer'='layer', label='extra-morphs'),
                           MorphInfo(name='name', mtype='mtype', 'layer'='layer', label='extra-morphs')])
    morph2 = all_morphs + more_morphs

Writing neurondb to disk.
~~~~~~~~~~~~~~~~~~~~~~~~~

The data can be written to disk with both the XML and DAT format.

.. code:: python

    MorphDB.from_neurondb('path1').write('neurondb.dat')

Analysing morph-stats features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NeuroM morphometrics can be extracted with ``self.features(config)``
where ``config`` is a `morph-stats
configuration <https://neurom.readthedocs.io/en/latest/morph_stats.html>`__.

Labels can be used to compare morphometrics across different datasets.

.. code:: python

    db = MorphDB.from_neurondb('path1', label='dataset-1')
    db += MorphDB.from_neurondb('path2', label='dataset-2')
    db += MorphDB.from_neurondb('path3', label='dataset-3')
    db += MorphDB.from_neurondb('path4', label='dataset-4')

    features = db.features(config)

    for dataset, df in features.groupby(('neuron', 'label')):
        print(dataset, df)
