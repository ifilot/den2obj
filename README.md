# Den2Obj

[![C/C++ CI](https://github.com/ifilot/den2obj/actions/workflows/build.yml/badge.svg)](https://github.com/ifilot/den2obj/actions/workflows/build.yml)
[![C/C++ CI](https://github.com/ifilot/den2obj/actions/workflows/build-openvdb.yml/badge.svg)](https://github.com/ifilot/den2obj/actions/workflows/build-openvdb.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Purpose
Converts VASP density files (i.e. CHGCAR / PARCHG) or a Gaussian cube file to a [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file), [Stanford .ply file](https://en.wikipedia.org/wiki/PLY_(file_format)) or [OpenVDB format](https://www.openvdb.org/).

## Example image
![3D Reaction-Diffusion system](img/reac_diff_3d_network_small.png)

*The isosurface above represents the concentration profile of a reaction-diffusion system in 3D using Gray-Scott kinetics. The isosurface has been generated using den2obj and rendered using [Blender](https://www.blender.org/).*

## Compilation instructions

### Debian Latest

Getting the dependencies
```
sudo apt install build-essential cmake libtclap-dev libboost-all-dev libopenvdb-dev libtbb-dev pkg-config libcppunit-dev libeigen3-dev
```

To compile, run the following commands:
```
git clone https://github.com/ifilot/den2obj.git
cd den2obj
mkdir build
cd build
cmake -DMOD_OPENVDB=1 ../src
make -j5
```

### Ubuntu Latest

The stable OpenVDB library (`libopenvdb`) under Ubuntu is incompatible with the Threading Building Blocks (`libtbb`) library. To solve this, manually compile and install OpenVDB 8.2 using the following instructions.

```
get https://github.com/AcademySoftwareFoundation/openvdb/archive/refs/tags/v8.2.0.tar.gz && tar -xvzf v8.2.0.tar.gz
mkdir openvdb-build && cd openvdb-build && cmake ../openvdb-8.2.0 -DCMAKE_INSTALL_PREFIX=/opt/openvdb && make -j9 && sudo make install
```

Thereafter, clone, configure and compile Den2Obj and link against OpenVDB 8.2.

```
git clone https://github.com/ifilot/den2obj.git
cd den2obj
mkdir build
cd build
cmake -DMOD_OPENVDB=1 ../src
make -j5
```

## Usage

### Isosurfaces

```
<path to>/den2obj -i CHGCAR -o <filename.obj> -v <isovalue>
```

Example:
```
./den2obj -i CHGCAR -o orbital.obj -v 0.1
```

## Options

* `-c`: Center the structure, i.e. the center of the structure is placed at the origin of the coordinate system.
* `-t`: Converts one file format to another. File formats are auto-recognized based on the extensions.

### Conversions

Converting CHGCAR to OpenVDB
```
./den2obj -i CHGCAR_xxx -o xxx.vdb -t
```

Converting CHGCAR to D2O
```
./den2obj -i CHGCAR_xxx -o xxx.d2o -t
```

Converting CHGCAR to D2O
```
./den2obj -i CHGCAR_xxx -o xxx.d2o -t
```

Supported input types:
* CHGCAR
* PARCHG
* LOCPOT
* Gaussian cube (.cub)
* D2O files (.d2o)

Supported dense output types:
* D2O
* OpenVDB

Supported isosurface object types:
* [Stanford .ply file](https://en.wikipedia.org/wiki/PLY_(file_format))
* [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file)

## D2O file format

The D2O file format is native to `Den2Obj`. This file format stores the scalarfield
in binary format and uses compression to generate small files which are fast to
read from. More information on the file format can be found in the
[documentation](https://den2obj.imc-tue.nl).
