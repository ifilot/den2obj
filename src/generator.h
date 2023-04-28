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

#ifndef _GENERATOR
#define _GENERATOR

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include "d2o_format.h"
#include "generator_benzene_data.h"
#include "generator_sto3g_data.h"

/**
 * Auxiliary class that can build ScalarFields
 **/
class Generator {
private:

    enum class GeneratorDataset {
        DS_GENUS2,
        DS_BENZENE_HOMO,
    };

    const std::unordered_map<std::string, GeneratorDataset> datasets = {
        {"genus2", GeneratorDataset::DS_GENUS2},
        {"benzene_homo", GeneratorDataset::DS_BENZENE_HOMO},
    };

public:
    Generator();

    void build_dataset(const std::string& name,
                       const std::string& filename) const;

private:
    std::vector<fpt> genus2(fpt sz, size_t npts) const;

    std::vector<fpt> benzene_molecular_orbital(fpt sz, size_t npts,
                                               unsigned int mo_id) const;

    fpt calculate_gto_normalization_constant(fpt alpha,
                                             unsigned int l,
                                             unsigned int m,
                                             unsigned int n) const;

    fpt calculate_gto_amp(const Vec3& pos, const Vec3& apos,
                          fpt alpha, unsigned int l,
                          unsigned int m, unsigned int n) const;

    fpt calculate_mo_amp(const Vec3& pos, unsigned int mo_id) const;
};

#endif // _GENERATOR
