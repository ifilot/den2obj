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
        DS_BENZENE_LUMO,
    };

    const std::unordered_map<std::string, GeneratorDataset> datasets = {
        {"genus2", GeneratorDataset::DS_GENUS2},
        {"benzene_homo", GeneratorDataset::DS_BENZENE_HOMO},
        {"benzene_lumo", GeneratorDataset::DS_BENZENE_LUMO},
    };

public:
    /**
     * @brief      Constructs a new instance.
     */
    Generator();

    /**
     * @brief      Build a dataset
     *
     * @param[in]  name      Name of the dataset
     * @param[in]  filename  File to output dataset to
     * @param[in]  algo_id   Compression algorithm
     */
    void build_dataset(const std::string& name,
                       const std::string& filename,
                       unsigned int algo_id = 0) const;

private:
    /**
     * @brief      Generator for the `genus2` dataset
     *
     * @param[in]  sz    Size in each cartesian direction
     * @param[in]  npts  Number of points in each direction
     *
     * @return     Vector containing datapoints for genus2 dataset
     */
    std::vector<fpt> genus2(fpt sz, size_t npts) const;

    /**
     * @brief      Generator for benzene molecular orbitals
     *
     * @param[in]  sz     Size in each cartesian direction
     * @param[in]  npts   Number of points in each direction
     * @param[in]  mo_id  Which molecular orbital (0-35) to produce
     *
     * @return     Vector containing datapoints for benzene molecule orbital
     */
    std::vector<fpt> benzene_molecular_orbital(fpt sz, size_t npts,
                                               unsigned int mo_id) const;

    /**
     * @brief      Calculate Gaussian Type Orbital Normalization Constant
     *
     * @param[in]  alpha  Exponent
     * @param[in]  l      Power in x-direction
     * @param[in]  m      Power in y-direction
     * @param[in]  n      Power in z-direction
     *
     * @return     The gto normalization constant.
     */
    fpt calculate_gto_normalization_constant(fpt alpha,
                                             unsigned int l,
                                             unsigned int m,
                                             unsigned int n) const;

    /**
     * @brief      Calculate the amplitude (value) of a GTO
     *
     * @param[in]  pos    Position to evaluate
     * @param[in]  apos   Location of the GTO
     * @param[in]  alpha  Exponent
     * @param[in]  l      Power in x-direction
     * @param[in]  m      Power in y-direction
     * @param[in]  n      Power in z-direction
     *
     * @return     The gto amp.
     */
    fpt calculate_gto_amp(const Vec3& pos, const Vec3& apos,
                          fpt alpha, unsigned int l,
                          unsigned int m, unsigned int n) const;

    /**
     * @brief      Calculate the ampltitude (value) of the molecular orbital
     *
     * @param[in]  pos    Position in three-dimensional space
     * @param[in]  mo_id  Which molecular orbital to produce
     *
     * @return     molecular orbital amplitude
     */
    fpt calculate_mo_amp(const Vec3& pos, unsigned int mo_id) const;
};

#endif // _GENERATOR
