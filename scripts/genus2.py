# -*- coding: utf-8 -*-

#
# Generate the "genus-2" scalar field. The function has
# been taken from: https://doi.org/10.1109/3DV.2015.36
#

import numpy as np
import gzip

def main():
    dim = 2
    grid = generate_grid(dim = dim)
    store_d2o(grid,dim)

def store_d2o(grid, dim):
    """
    Store Genus2 scalar field as a d2o file. See
    https://den2obj.imc-tue.nl/d2o_fileformat.html
    for the specifications of the d2o file type.
    """
    f = open('genus2.d2o', 'wb')  
    
    # write header
    f.write('D2O'.encode('utf-8'))
    
    # write protocol identifier
    f.write(np.uint32(1).tobytes())
    
    # write unit cell matrix
    unitcell = np.diag(np.ones(3, dtype=np.float32) * dim * 2)
    f.write(unitcell.flatten().tobytes())
    
    # write grid dimensions
    f.write(np.array(grid.shape, dtype=np.uint32).tobytes())
    
    # write floating point bytesize
    f.write(np.uint8(4).tobytes())
    
    # write data stream
    compressed_data = gzip.compress(np.array(grid.flatten(), dtype=np.float32).tobytes(),
                                    compresslevel=9)
    f.write(np.uint64(len(compressed_data)).tobytes())
    f.write(compressed_data)
    
    f.close()
    
def generate_grid(dim = 2, npts = 100):
    x = np.linspace(-dim, dim, npts)
    y = np.linspace(-dim, dim, npts)
    z = np.linspace(-dim, dim, npts)
    
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    
    return genus2(xx,yy,zz)
    
def genus2(x,y,z):
    F = 2*y*(y**2 - 3*x**2) * (1 - z**2) + \
        (x**2 + y**2)**2 - (9*z**2 - 1) * (1 - z**2)
    return F
    
if __name__ == '__main__':
    main()