# Den2Obj

## Purpose
Converts VASP density files (i.e. CHGCAR / PARCHG) or a (custom-format) [binary file](#binary-source) to a [Wavefront .obj file](https://en.wikipedia.org/wiki/Wavefront_.obj_file) or a [Stanford .ply file](https://en.wikipedia.org/wiki/PLY_(file_format)).

## Example image
![3D Reaction-Diffusion system](img/reac_diff_3d_network_small.png)

*The isosurface above represents the concentration profile of a reaction-diffusion system in 3D using Gray-Scott kinetics. The isosurface has been generated using den2obj and rendered using [Blender](https://www.blender.org/).*

## Compilation instructions

Getting the dependencies
```
sudo apt install build-essential cmake libglm-dev libtclap-dev libboost-all-dev
```

To compile, run the following commands:
```
git clone https://github.com/ifilot/den2obj.git
cd den2obj
mkdir build
cd build
cmake ../src
make -j5
```

## Usage

```
<path to>/den2obj -i CHGCAR -o <filename.obj> -v <isovalue>
```

Example:
```
./den2obj -i CHGCAR -o orbital.obj -v 0.1
```

Example output:
```
--------------------------------------------------------------
Executing DEN2OBJ v.0.4.0
Author: Ivo Filot <i.a.w.filot@tue.nl>
Website: https://github.com/ifilot/den2obj
--------------------------------------------------------------
Using isovalue: 0.1
Identified 1608 faces.
Writing to orbital.obj
--------------------------------------------------------------
Done in 0.0177371 seconds.
```

## Options

* `-c`: Center the structure, i.e. the center of the structure is placed at the origin of the coordinate system.
* `-p`: Write output as a binary `.ply` file rather than `.obj` file. The program automatically detects the [endianness](https://en.wikipedia.org/wiki/Endianness) of your system.
* `-b`: Read from binary source rather than `CHGCAR` or `PARCHG` file (read about the format of this [binary file](#binary-source) below!).

## Binary source
Instead of a `CHGCAR` or `PARCHG` file, a binary file can be supplied which should be formatted by first stating the number of grid points in x, y and z direction and the size of the floating point variable (only float and double are currently supported) followed by a list of 64 bit floats for the value at each grid point. Herein, it assumed that "z" is the slowest index and "x" the fastest index.

In short:
```
<size x><size y><size z><size float><...data...>
```
