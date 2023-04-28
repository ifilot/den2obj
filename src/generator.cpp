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

void Generator::build_dataset(const std::string& name, const std::string& filename) const {
    auto got = this->datasets.find(name);
    if(got == this->datasets.end()) {
        throw std::runtime_error("Invalid dataset: " + name);
    } else {
        switch(got->second) {
            case GeneratorDataset::DS_GENUS2:
            {
                const unsigned int npts = 100;
                const fpt sz = 2.0f;
                std::array<unsigned int, 3> grid_dimensions = {npts, npts, npts};
                std::vector<fpt> grid = this->genus2(sz, npts);
                MatrixUnitcell mat = MatrixUnitcell::Zero();
                for(unsigned int i=0; i<3; i++) {
                    mat(i,i) = 2.0 * sz;
                }
                D2OFormat::write_d2o_file(filename, grid, grid_dimensions, mat);
            }
            break;
            case GeneratorDataset::DS_BENZENE_HOMO:
            {
                const unsigned int npts = 150;
                const fpt sz = 6.0f;
                std::array<unsigned int, 3> grid_dimensions = {npts, npts, npts};
                std::vector<fpt> grid = this->benzene_molecular_orbital(sz, npts, 20);
                MatrixUnitcell mat = MatrixUnitcell::Zero();
                for(unsigned int i=0; i<3; i++) {
                    mat(i,i) = 2.0 * sz;
                }
                D2OFormat::write_d2o_file(filename, grid, grid_dimensions, mat);
            }
            break;
            default:
                throw std::logic_error("Cannot execute dataset generation.");
            break;
        }
    }
}

std::vector<fpt> Generator::genus2(fpt sz, size_t npts) const {
    std::vector<fpt> xx(npts);

    // build grid points in one nptsensions
    for(unsigned int i=0; i<npts; i++) {
        xx[i] = -sz + (2 * sz / (fpt)(npts)) * i;
    }

    // build grid
    std::vector<fpt> f(npts*npts*npts);
    #pragma omp parallel for schedule(static)
    for(unsigned int i=0; i<npts; i++) { // loop over z
        const float z = xx[i];
        for(unsigned int j=0; j<npts; j++) { // loop over y
            const float y = xx[j];
            for(unsigned int k=0; k<npts; k++) { // loop over x
                const float x = xx[k];
                const unsigned int idx = i * npts * npts + j * npts + k;
                f[idx] = 2 * y * (y*y - 3 * x*x) * (1 - z*z) +
                         (x*x + y*y) * (x*x + y*y) - (9*z*z - 1) * (1 - z*z);
            }
        }
    }

    return f;
}

std::vector<fpt> Generator::benzene_molecular_orbital(fpt sz, size_t npts,
                                                      unsigned int mo_id) const {
    std::vector<fpt> xx(npts);

    // build grid points in one nptsensions
    for(unsigned int i=0; i<npts; i++) {
        xx[i] = -sz + (2 * sz / (fpt)(npts)) * i;
    }

    // build grid
    std::vector<fpt> f(npts*npts*npts);
    #pragma omp parallel for schedule(static)
    for(unsigned int i=0; i<npts; i++) { // loop over z
        const float z = xx[i];
        for(unsigned int j=0; j<npts; j++) { // loop over y
            const float y = xx[j];
            for(unsigned int k=0; k<npts; k++) { // loop over x
                const float x = xx[k];
                const unsigned int idx = i * npts * npts + j * npts + k;
                f[idx] = this->calculate_mo_amp(Vec3(x,y,z), mo_id);
            }
        }
    }

    return f;
}

fpt Generator::calculate_gto_normalization_constant(fpt alpha,
                                                    unsigned int l,
                                                    unsigned int m,
                                                    unsigned int n) const {

    static const fpt pi = boost::math::constants::pi<fpt>();

    fpt nom =    std::pow(2.0, 2.0 * (l + m + n) + 3.0 / 2.0) *
                 std::pow(alpha, (l + m + n) + 3.0 / 2.0);

    fpt denom = (l < 1 ? 1 : boost::math::double_factorial<fpt>(2 * l - 1) ) *
                (m < 1 ? 1 : boost::math::double_factorial<fpt>(2 * m - 1) ) *
                (n < 1 ? 1 : boost::math::double_factorial<fpt>(2 * n - 1) ) *
                std::pow(pi, 3.0 / 2.0);

    return std::sqrt(nom / denom);
}

fpt Generator::calculate_gto_amp(const Vec3& pos, const Vec3& apos,
                                 fpt alpha, unsigned int l,
                                 unsigned int m, unsigned int n) const {

    fpt r2 = (pos - apos).squaredNorm();

    const fpt norm = this->calculate_gto_normalization_constant(alpha, l, m, n);

    return norm * std::pow(pos[0] - apos[0], l) *
                  std::pow(pos[1] - apos[1], m) *
                  std::pow(pos[2] - apos[2], n) *
                  std::exp(-alpha * r2);
}

fpt Generator::calculate_mo_amp(const Vec3& pos, unsigned int mo_id) const {
    fpt val = 0.0;

    for(unsigned int i=0; i<6; i++) { // loop over C atoms
        const Vec3 apos(GeneratorData::benzene_atompos[i][0] * 1.88973f,
                        GeneratorData::benzene_atompos[i][1] * 1.88973f,
                        GeneratorData::benzene_atompos[i][2] * 1.88973f);

        for(unsigned int j=0; j<5; j++) { // loop over C basis functions
            const unsigned int cgf_id = std::min((unsigned int)2,j);
            const float lcao_coeff = GeneratorData::benzene_orbc[5 * i + j][mo_id];
            unsigned int l = 0;
            unsigned int m = 0;
            unsigned int n = 0;

            if(j == 2) {
                l = 1;
            } else if(j == 3) {
                m = 1;
            } else if(j == 4) {
                n = 1;
            }

            for(unsigned int k=0; k<3; k++) { // loop over GTOs
                const float alpha = GeneratorData::basis_c[cgf_id][k][0];
                const float gto_coeff = GeneratorData::basis_c[cgf_id][k][1];
                const float amp = calculate_gto_amp(pos, apos, alpha, l, m, n);

                val += lcao_coeff * gto_coeff * amp;
            }
        }
    }

    for(unsigned int i=0; i<6; i++) { // loop over H atoms
        const Vec3 apos(GeneratorData::benzene_atompos[6 + i][0] * 1.88973f,
                        GeneratorData::benzene_atompos[6 + i][1] * 1.88973f,
                        GeneratorData::benzene_atompos[6 + i][2] * 1.88973f);

        const float lcao_coeff = GeneratorData::benzene_orbc[30 + i][mo_id];
        unsigned int l = 0;
        unsigned int m = 0;
        unsigned int n = 0;

        for(unsigned int k=0; k<3; k++) { // loop over GTOs
            const float alpha = GeneratorData::basis_h[k][0];
            const float gto_coeff = GeneratorData::basis_h[k][1];
            const float amp = calculate_gto_amp(pos, apos, alpha, l, m, n);

            val += lcao_coeff * gto_coeff * amp;
        }
    }

    return val;
}
