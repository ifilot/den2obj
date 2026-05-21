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
* `Zstandard <https://facebook.github.io/zstd/>`_ (zstd data compression)
* `Blosc <https://www.blosc.org/>`_ (blocked and shuffled data compression)
* `CPPUnit <https://sourceforge.net/projects/cppunit/>`_ (unit testing)

.. note::
   * The instructions covered in this guide assume you are running a  
     Debian-based Linux distro such as Debian, Ubuntu, or Kali. 
   * If you are running Windows and would like to use :program:`Den2Obj`, one 
     option is to use `Debian for Windows Subsystem for Linux (WSL) <https://apps.microsoft.com/store/detail/debian/9MSVKQC78PK6>`_.
     The compilation instructions below can be readily used.

To ensure that all the packages are installed, one can run the following::

    sudo apt install build-essential cmake libtclap-dev libboost-all-dev \
    pkg-config libcppunit-dev libeigen3-dev liblzma-dev zlib1g-dev libbz2-dev \
    libzstd-dev libblosc-dev

Standard compilation
--------------------

Compilation of :program:`Den2Obj` is fairly straightforward and a typical procedure
looks as follows::

    git clone https://github.com/ifilot/den2obj.git
    cd den2obj
    mkdir build && cd build
    cmake ../src
    cmake --build . --parallel

To install :program:`Den2Obj`, you can in addition run::

    sudo cmake --install .

which will place a copy of the ``den2obj`` executable in ``/usr/local/bin/den2obj``.

Testing
-------

To test :program:`Den2Obj`, one can run the following after compilation::

	ctest --output-on-failure

A successful test run should report that all tests passed. The exact number of
tests and their timings can differ between standard and coverage builds.

If the test is for some reason failing, one can run the following to produce
more output::

    ctest --output-on-failure

.. note::

    If the tests are continously failing for you, you are warmly invited
    to `open an issue on the Github page <https://github.com/ifilot/den2obj/issues>`_.
