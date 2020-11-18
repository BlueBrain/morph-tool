MorphTool
=========

A toolbox for morphology editing. It aims to provide small helper programs that perform simple tasks.

Currently MorphTool provides:

- A morphology diffing tool (via CLI or python)
- A file converter: to convert morphology files from/to the following formats: SWC, ASC, H5
- A soma area calculator (as computed by NEURON, requires the NEURON python module)

Installation
------------

In a fresh virtualenv:

.. code:: bash

    pip install morph-tool

Usage
-----

In a shell, do:

.. code:: bash

    morph-tool --help

Currently the three sub-command are:

.. code:: bash

    morph-tool convert input_file output_file
    morph-tool diff morph1 morph2
    morph-tool soma-surface input_file

Morphology diffing
------------------
One can compare two morphologies with the CLI

.. code:: bash

   morph-tool diff morph1 morph2

whose error code is 0 if morphologies are the same, else 1.

Morphologies with different formats can be compared.

Morphologies are considered different if one or more of the following properties differ:

- number of root sections
- sections type
- sections point array
- sections diameter array
- sections perimeter array
- sections number of children

The soma are *not* taken into consideration

The same functionality is also avalaible through the python API:

.. code:: python

   from morph_tool import diff

   # The result can be used as a boolean:
   if diff(filename1, filename2):
       print('morphologies differ')

   # And also contains information about how morphologies differ
   result = diff(filename1, filename2)
   print(result.info)



Converter
---------

What is supported ?
~~~~~~~~~~~~~~~~~~~

The converter can be used to write morphology in different formats.
Currently, supported formats are ASC, SWC and H5.

As each format has its own specificities, data specific to a given
format will be discarded. This means the following will be lost during
conversion:

* spines (present in the ASC format)

* all H5 metadata

* the perimeter and mitochondrial data of the `H5 format <https://bbpteam.epfl.ch/documentation/Morphology%20Documentation-0.0.2/h5v1.html>`__

Soma intricacies
~~~~~~~~~~~~~~~~

Multiple formats are being used to represent somas (mainly) depending on
the file format. For more information about file format, see the `neuromorpho.org specification <http://neuromorpho.org/SomaFormat.html>`__ or `MorphIO
specification <https://github.com/BlueBrain/MorphIO/blob/master/doc/specification.md#soma-formats>`__

Because different soma format represent soma in different planes, soma
format conversion is not a bijective transformation. For example, it is
not possible to have an accurate conversion from a soma contour in the
XY plane from a H5 file to a SWC soma which is represented as a cylinder
along Y.

As a result, it has been chosen to take the soma surface as an
invariant. This means soma surfaces of the input and output morphologies, as computed by NEURON, should be preserved.

Here are the possible cases for the soma conversion:

- SWC input file:

  - SWC output file -> no conversion
  - H5 or ASC output file:

    Depending on the original soma type:

    - Soma stack of cylinders:
      The soma is converted to a contour in the XY plane.
      The points of the new contour are the outline of the soma stack projected in the XY plane.
    - Soma three point cylinder:
      The soma becomes a sphere of same surface. The contour made by the circle of biggest section in the XY plane is sampled in 20 points written to disk.
    - Soma sphere (soma represented by a single point representing the center of a sphere and its radius): the contour made by the circle of biggest section in the XY plane is sampled in 20 points written to disk.
    - other:
      Not in SWC spec -> not supported

- H5 or ASC input file:

  - H5 output file -> no conversion needed
  - ASC output file.

    Depending on soma type:

    - Soma single point sphere (soma represented by a single point representing the center of a sphere and its radius): the contour made by the circle of biggest section in the XY plane is sampled in 20 points written to disk.
    - Soma contour: no conversion needed
    - other: not in H5/ASC specs -> not supported

  - SWC:

    Depending on soma format:

    - Soma single point sphere: no conversion needed
    - Soma contour: A soma stack of cylinder is generated.
      Each cylinder of the stack has its center and its axis along the principal direction of the contour.
      The radius of each stack is chosen such that it minimises the distance between the cylinder and the contour.
    - other: not in H5/ASC specs -> not supported

Example:

.. code:: python

   from morph_tool import convert
   convert(inputfile, outputfile)

   # Additionally the morphology can be recentered or written according to the NEURON neurite order during the conversion
   convert(inputfile, outputfile, recenter=True, nrn_order=True)

The same for bash

.. code:: bash

   morph-tool convert file ./inputfile ./outputfile
   # with additional options
   morph-tool convert file --recenter --nrn-order ./inputfile ./outputfile
   # or an entire folder
   morph-tool convert folder -ext SWC ./h5_input_folder ./swc_output_folder
   # for more info use
   morph-tool convert folder --help

NRN simulator compartment coordinates
-------------------------------------

The NRN simulator splits each section into chunks of equal length (equal only among a given section).
These compartments do not really exist in the physical world but we can remap them to paths
along the section. Each compartment can be associated to a path (a list of 3D points) such
that the path and the compartment have the same path-length.

The following function can be used to access the mapping NeuroM section ID -> list of paths for the section:

.. code:: python

          morph_tool.nrnhines.NeuroM_section_to_NRN_compartment_paths


Example (in 2D) for one section:

.. code::

                   (1, 2) ------ (2, 2)
                      |
                      |
                      |
                      |
                      |
                      |
                      |
                      |
                      |
    (0, 0) ------- (1, 0)


Splitting this section into 3 compartments would results in the following paths:

1.

.. code::

    [[0.        , 0.        ],
     [1.        , 0.        ],
     [1.        , 0.33333333]]

2.

.. code::

   [[1.        , 0.33333333],
    [1.        , 1.66666667]]

3.

.. code::

   [[1.        , 1.66666667],
    [1.        , 2.        ],
    [2.        , 2.        ]]


Dendrogram with synapses
------------------------

This functionality is available only when the package is installed with **dendrogram** extras:

.. code:: bash

    pip install morph-tool

Draw NeuroM dendrogram with synapses on it. Synapses must be represented as a DataFrame. Required
columns in this dataframe are:

.. code:: python

    from morph_tool import dendrogram
    required_columns = [dendrogram.SOURCE_NODE_ID, dendrogram.TARGET_NODE_ID,
                        dendrogram.POST_SECTION_ID, dendrogram.POST_SECTION_POS,
                        dendrogram.PRE_SECTION_ID, dendrogram.PRE_SECTION_POS]

or equivalently

.. code:: python

    required_columns = ['@source_node', '@target_node',
                        'afferent_section_id', 'afferent_section_pos',
                        'efferent_section_id', 'efferent_section_pos']


For usage examples look at ``examples/dendrogram.py``.

Contributing
------------

If you want to improve the project or you see any issue, every contribution is welcome.
Please check the `contribution guidelines <CONTRIBUTING.md>`__ for more information.

Acknowledgements
----------------

This research was supported by the EBRAINS research infrastructure, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under the Specific Grant Agreement No. 945539 (Human Brain Project SGA3).

License
-------

morph-tool is licensed under the terms of the GNU Lesser General Public License version 3.
Refer to COPYING.LESSER and COPYING for details.
