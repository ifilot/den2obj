/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   EDP is free software:                                                *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   EDP is distributed in the hope that it will be useful,               *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _D2O_FORMAT
#define _D2O_FORMAT

#include <vector>
#include <iostream>
#include <array>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include "math.h"

/**
 * @brief      Write a D2O file
 *
 * @param[in]  filename         The filename
 * @param[in]  protocol_id      The protocol identifier
 * @param[in]  gridptr          Vector containing the grid points
 * @param      grid_dimensions  Dimensions in each unit cell direction
 * @param[in]  mat              Unitcell matrix
 */
void write_d2o_file(const std::string& filename,
                    uint32_t protocol_id,
                    const std::vector<fpt>& gridptr,
                    std::array<unsigned int, 3>& grid_dimensions,
                    const MatrixUnitcell& mat);

/**
 * @brief      Compress data stream using all possible compression algos
 *
 * @param[in]  origin  Uncompressed data stream
 *
 * @return     Vector containing compressed streams as strings
 */
std::vector<std::string> d2o_compress_all(const std::istringstream& origin);

#endif // _D2O_FORMAT