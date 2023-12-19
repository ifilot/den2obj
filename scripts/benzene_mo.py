from pyqint import PyQInt, cgf, gto, Molecule, HF
from pytessel import PyTessel
import numpy as np

mol = Molecule()
mol.add_atom('C',  0.0000000015, -1.3868467444, 0.0000000000, unit='angstrom')
mol.add_atom('C',  1.2010445126, -0.6934233709, 0.0000000000, unit='angstrom')
mol.add_atom('C',  1.2010445111,  0.6934233735, 0.0000000000, unit='angstrom')
mol.add_atom('C', -0.0000000015,  1.3868467444, 0.0000000000, unit='angstrom')
mol.add_atom('C', -1.2010445126,  0.6934233709, 0.0000000000, unit='angstrom')
mol.add_atom('C', -1.2010445111, -0.6934233735, 0.0000000000, unit='angstrom')
mol.add_atom('H',  0.0000000027, -2.4694205285, 0.0000000000, unit='angstrom')
mol.add_atom('H',  2.1385809117, -1.2347102619, 0.0000000000, unit='angstrom')
mol.add_atom('H',  2.1385809090,  1.2347102666, 0.0000000000, unit='angstrom')
mol.add_atom('H', -0.0000000027,  2.4694205285, 0.0000000000, unit='angstrom')
mol.add_atom('H', -2.1385809117,  1.2347102619, 0.0000000000, unit='angstrom')
mol.add_atom('H', -2.1385809090, -1.2347102666, 0.0000000000, unit='angstrom')

results = HF().rhf(mol, 'sto3g')

for row in results['orbc']:
    print("\t{", end='')
    for val in row:
        print("%+12.6f, " % val, end='')
    print("},")
    
def build_isosurface(filename, cgfs, coeff, isovalue):
    # generate some data
    dz = 10.0
    sz = 150
    integrator = PyQInt()
    grid = integrator.build_rectgrid3d(-dz, dz, sz)
    scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (sz, sz, sz))
    unitcell = np.diag(np.ones(3) * dz * 2)

    pytessel = PyTessel()
    vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
    pytessel.write_ply(filename, vertices, normals, indices)
    
build_isosurface('benzene_p.ply', results['cgfs'], results['orbc'][:,20], 0.01)
build_isosurface('benzene_n.ply', results['cgfs'], results['orbc'][:,20], -0.01)