.. _installation:
.. index:: Installation

Installation
============

.. contents::
   :local:

Prerequisites
-------------

:program:`Den2Obj` is developed for Linux operating systems. In order to 
compile :program:`Den2Obj` on your system, you need to ensure the following 
libraries are available to you:

* `Eigen3 <https://eigen.tuxfamily.org>`_ (matrix algebra)
* `Boost <https://www.boost.org/>`_ (common routines)
* `TCLAP <https://tclap.sourceforge.net/>`_ (command line instruction library)
* `BZIP2 <https://sourceware.org/bzip2/>`_ (bzip2 data compression)
* `GZIP <https://www.gnu.org/software/gzip/>`_ (gzip data compression)
* `LZMA <https://7-zip.org/>`_ (lzma data compression)
* `CPPUnit <https://sourceforge.net/projects/cppunit/>`_ (unit testing)

.. note::
   * The instructions covered in this guide assume you are running a  
     Debian-based Linux distro such as Debian, Ubuntu, or Kali. 
   * If you are running Windows and would like to use :program:`Den2Obj`, one 
     option is to use `Debian for Windows Subsystem for Linux (WSL) <https://apps.microsoft.com/store/detail/debian/9MSVKQC78PK6>`_.
     The compilation instructions below can be readily used.

To ensure that all the packages are installed, one can run the following::

    sudo apt install build-essential cmake libtclap-dev libboost-all-dev \ 
    pkg-config libcppunit-dev libeigen3-dev liblzma-dev zlib1g-dev libbz2-dev

Standard compilation
--------------------

Compilation of :program:`Den2Obj` is fairly straightforward and a typical procedure
looks as follows::

    git clone https://github.com/ifilot/den2obj.git
    cd den2obj
    mkdir build && cd build
    cmake ../src
    make -j9

To install :program:`Den2Obj`, you can in addition run::

    sudo make install

which will place a copy of the ``den2obj`` executable in ``/usr/local/bin/den2obj``.

Compilation with OpenVDB module
-------------------------------

Using the `OpenVDB <https://www.openvdb.org/>`_ library, it is possible
to convert a density object to an OpenVDB file which can be used to render
the density clouds in `Blender <https://www.blender.org/>`_. Compilation
of :program:`Den2Obj` using this module requires supplying an additional
argument to CMake.

.. warning::
    The OpenVDB package in the latest Ubuntu LTS (22.04 LTS) is incompatible with the Thread Building Blocks (tbb) library. As a workaround, one is
    required to compile a newer OpenVDB. OpenVDB v8.2 has been tried and
    tested by us. If you are compiling for Ubuntu 22.04 LTS, please read the
    additional instructions to compile OpenVDB.

Ensure that the following libraries are installed by running::

	sudo apt install build-essential cmake libtclap-dev libboost-all-dev \
	libopenvdb-dev libtbb-dev pkg-config libcppunit-dev libeigen3-dev \
    liblzma-dev zlib1g-dev libbz2-dev

Compilation :program:`Den2Obj` with the OpenVDB module is as follows::

    git clone https://github.com/ifilot/den2obj.git
    cd den2obj
    mkdir build && cd build
    cmake -DMOD_OPENVDB=1 ../src
    make -j9

Testing
-------

To test :program:`Den2Obj`, one can run the following after compilation::

	make test

A succesfull test should produce an output similar to the one found below::

    Running tests...
    Test project /mnt/c/PROGRAMMING/CPP/den2obj/build
        Start 1: DatasetSetup
    1/6 Test #1: DatasetSetup .....................   Passed    2.49 sec
        Start 3: TestIsosurface
    2/6 Test #3: TestIsosurface ...................   Passed    1.07 sec
        Start 4: TestScalarField
    3/6 Test #4: TestScalarField ..................   Passed    0.39 sec
        Start 5: TestD2OFileFormat
    4/6 Test #5: TestD2OFileFormat ................   Passed    0.02 sec
        Start 6: TestGenerator
    5/6 Test #6: TestGenerator ....................   Passed    8.34 sec
        Start 2: DatasetCleanup
    6/6 Test #2: DatasetCleanup ...................   Passed    0.00 sec

    100% tests passed, 0 tests failed out of 6

    Total Test time (real) =  12.45 sec

If the test is for some reason failing, one can run the following to produce
more output::

    CTEST_OUTPUT_ON_FAILURE=TRUE make test

.. note::

    If the tests are continously failing for you, you are warmly invited
    to `open an issue on the Github page <https://github.com/ifilot/den2obj/issues>`_.
