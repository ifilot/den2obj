/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   DEN2OBJ is free software:                                            *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DEN2OBJ is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include "generator.h"

Generator::Generator() {}

void Generator::build_dataset(const std::string& name, const std::string& filename) {
    auto got = this->datasets.find(name);
    if(got == this->datasets.end()) {
        throw std::runtime_error("Invalid dataset: " + name);
    } else {
        switch(got->second) {
            case GeneratorDataset::DS_GENUS2:
            {
                const unsigned int dim = 100;
                const fpt sz = 2.0f;
                std::array<unsigned int, 3> grid_dimensions = {dim, dim, dim};
                std::vector<fpt> grid = this->genus2(sz, dim);
                MatrixUnitcell mat = MatrixUnitcell::Zero();
                for(unsigned int i=0; i<3; i++) {
                    mat(i,i) = 2.0 * sz;
                }
                write_d2o_file(filename, 2, grid, grid_dimensions, mat);
            }
            break;
            default:
                throw std::logic_error("Cannot execute dataset generation.");
            break;
        }
    }
}

std::vector<fpt> Generator::genus2(fpt sz, size_t dim) {
    std::vector<fpt> xx(dim);

    // build grid points in one dimensions
    for(unsigned int i=0; i<dim; i++) {
        xx[i] = -sz + (2 * sz / (fpt)(dim-1)) * i;
    }

    // build grid
    std::vector<fpt> f(dim*dim*dim);
    #pragma omp parallel for schedule(static)
    for(unsigned int i=0; i<dim; i++) { // loop over z
        const float z = xx[i];
        for(unsigned int j=0; j<dim; j++) { // loop over y
            const float y = xx[j];
            for(unsigned int k=0; k<dim; k++) { // loop over x
                const float x = xx[k];
                const unsigned int idx = i * dim * dim + j * dim + k;
                f[idx] = 2 * y * (y*y - 3 * x*x) * (1 - z*z) +
                         (x*x + y*y) * (x*x + y*y) - (9*z*z - 1) * (1 - z*z);
            }
        }
    }

    return f;
}
