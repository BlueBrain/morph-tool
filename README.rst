MorphTool
=========

A toolbox for morphology editing. It aims to provide small helper programs that perform simple tasks.

Currently MorphTool provides:

- A morphology diffing tool (via CLI or python)
- A file converter: to convert morphology files from/to the following formats: SWC, ASC, H5
- A soma area calculator (as computed by NEURON, requires the NEURON python module)

Installation
============

In a fresh virtualenv:

.. code:: bash

    pip install --index-url https://bbpteam.epfl.ch/repository/devpi/bbprelman/dev/+simple/ morph-tool

Usage
=====

In a shell, do:

.. code:: bash

    morph-tool --help

Currently the three sub-command are:

.. code:: bash

    morph-tool convert input_file output_file
    morph-tool diff morph1 morph2
    morph-tool soma-surface input_file

Morphology diffing
==================
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
   result = diff(filename1, filename2)

Converter
=========

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

As a result, it has been choosen to take the soma surface as an
invariant. This means soma surfaces of the input and output morphologies, as computed by NEURON, should be the preserved.

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
    - Soma sphere (soma represented by a single point represeting the center of a sphere and its radius): the contour made by the circle of biggest section in the XY plane is sampled in 20 points written to disk.
    - other:
      Not in SWC spec -> not supported

- H5 or ASC input file:

  - H5 output file -> no conversion needed
  - ASC output file.

    Depending on soma type:

    - Soma single point sphere (soma represented by a single point represeting the center of a sphere and its radius): the contour made by the circle of biggest section in the XY plane is sampled in 20 points written to disk.
    - Soma contour: no conversion needed
    - other: not in H5/ASC specs -> not supported

  - SWC:

    Depending on soma format:

    - Soma single point sphere: no conversion needed
    - Soma contour: A soma stack of cylinder is generated.
      Each cylinder of the stack has its center and its axis along the principal direction of the contour.
      The radius of each stack is choosen such that it minimises the distance between the cylinder and the contour.
    - other: not in H5/ASC specs -> not supported
