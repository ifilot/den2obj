.. Den2Obj documentation master file, created by
   sphinx-quickstart on Sun Apr 16 08:39:47 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Den2Obj - A command line program for producing isosurfaces from density files
=============================================================================

.. image:: https://img.shields.io/github/v/tag/ifilot/den2obj?label=version
   :alt: GitHub tag (latest SemVer)
.. image:: https://github.com/ifilot/den2obj/actions/workflows/build.yml/badge.svg
   :target: https://github.com/ifilot/den2obj/actions/workflows/build.yml
.. image:: https://github.com/ifilot/den2obj/actions/workflows/build-openvdb.yml/badge.svg
   :target: https://github.com/ifilot/den2obj/actions/workflows/build-openvdb.yml
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

:program:`Den2Obj` is a command-line tool that construct isosurfaces from
densely packed scalar fields. :program:`Den2Obj` supports VASP charge files
such as CHGCAR and PARCHG, Gaussian .cube files as well as its own
:ref:`.d2o file format<D2O file format>`.

:program:`Den2Obj` can be used together with popular 3D rendering programs
such as `Blender <https://www.blender.org/>`_ to produce stunning visuals
of your research. Below, a few examples are shown. See also the
:ref:`Examples<Examples>` for a nice overview.

.. figure:: _static/img/ch4_valence_orbitals.png
   :alt: Canonical valence orbitals of CH4

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   background
   tutorial
   examples
   user_interface
   d2o_fileformat
   community_guidelines
   faq

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
