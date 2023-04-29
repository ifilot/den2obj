.. _userinterface:
.. index:: User Interface

User Interface
**************

:program:`Den2Obj` has two operational modes

1. Isosurface generation
2. Filetype conversion

Isosurface generation is the default mode and is used to generate isosurfaces
from scalar fields. Filetype conversion is mainly used to convert relatively
large storage formats such as CHGCAR or Gaussian cube files to the compressed
:ref:`.d2o file format<D2O file format>`.

Isosurface generation
=====================

To perform an isosurface creation, one simply runs::

    ./den2obj -i <path-to-scalarfield> -o <mesh-file> -v <isovalue> [-c]

wherein ``path-to-scalarfield`` is a scalar field file format of either
of the following types

* `CHGCAR <https://www.vasp.at/wiki/index.php/CHGCAR>`_, `PARCHG <https://www.vasp.at/wiki/index.php/PARCHG>`_ or `LOCPOT <https://www.vasp.at/wiki/index.php/LOCPOT>`_
* `Gaussian Cube file <https://gaussian.com/cubegen/>`_
* :ref:`Den2Obj .d2o file format<D2O file format>`

``mesh-file`` is any of the following supported formats

* `Stereolitography (.stl) file <https://en.wikipedia.org/wiki/STL_(file_format)>`_
* `Stanford (.ply) file <https://en.wikipedia.org/wiki/PLY_(file_format)>`_
* `Wavefront (.obj) file <https://en.wikipedia.org/wiki/Wavefront_.obj_file>`_

and ``isovalue`` the isovalue for the isosurface.

.. tip::

    For visualizing wave functions, two isosurfaces are needed. One for the
    positive part and one for the negative part of the wave function. Simply
    run :program:`Den2Obj` twice and put a minus sign in front of the isovalue.
    Take care to adjust the output file, else the original output file will
    be overwritten.

The argument ``-c`` is optional. If used, the isosurface will be centered at
the origin.

Filetype conversion
===================

:program:`Den2Obj` offers the conversion to two different file types.

* :ref:`Den2Obj .d2o file format<D2O file format>`
* `OpenVDB <https://www.openvdb.org/>`_

To perform the conversion, one executes::

    ./den2obj -i <path-to-scalarfield> -o <mesh-file> -t

.. note::

    Creating `OpenVDB` files requires the OpenVDB module to be compiled. See
    the :ref:`installation <Installation>` section for more details.

Dataset generation
==================

For testing and learning purposes, :program:`Den2Obj` can generate a number of datasets. This
is done via the ``-g <dataset name>`` directive, which takes a valid ``dataset name`` as input. 

For example, to build the ``genus2`` dataset, one runs::

    ./den2obj -g genus2 -o genus2.d2o

It is possible to overrule the preferred compression algorithm using the ``-a <compression algo>``
directive. For example, to force BZIP2 type of compression, one runs::

    ./den2obj -g genus2 -o genus2.d2o -a bzip2