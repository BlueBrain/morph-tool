MorphTool
=========

A toolbox for morphology editing. It aims at providing small programs
that perform simple tasks.

Currently MorphTool provides: - A file converter: to convert morphology
files from/to the following formats: SWC, ASC, H5

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

to list all functionalities, or:

.. code:: bash

    morph-tool convert input_file output_file

to convert a file.

Specifications
==============

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

As a result, it has been choosen to take the soma surface as an
invariant.

.. note:: This means soma surfaces of the input and output morphologies, as computed by NRN, should be the preserved.
