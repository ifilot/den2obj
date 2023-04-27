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

#include "d2o_format.h"

/**
 * Auxiliary class that can build ScalarFields
 **/
class Generator {
private:

    enum class GeneratorDataset {
        DS_GENUS2,
    };

    const std::unordered_map<std::string, GeneratorDataset> datasets = {
        {"genus2", GeneratorDataset::DS_GENUS2},
    };

public:
    Generator();

    void build_dataset(const std::string& name,
                       const std::string& filename);

private:
    std::vector<fpt> genus2(fpt sz, size_t dim);

};

#endif // _GENERATOR
