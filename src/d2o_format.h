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

#include <unordered_map>
#include <vector>
#include <iostream>
#include <array>
#include <fstream>

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

namespace D2OFormat {

    enum class CompressionAlgo {
        AUTO,
        GZIP,
        LZMA,
        BZIP2,
    };

    // list of compression algos
    static const std::unordered_map<std::string, CompressionAlgo> algos = {
        {"auto", CompressionAlgo::AUTO},
        {"gzip", CompressionAlgo::GZIP},
        {"lzma", CompressionAlgo::LZMA},
        {"bzip2", CompressionAlgo::BZIP2}
    };

    /**
     * @brief      Write a D2O file
     *
     * @param[in]  filename         The filename
     * @param[in]  gridptr          Vector containing the grid points
     * @param      grid_dimensions  Dimensions in each unit cell direction
     * @param[in]  mat              Unitcell matrix
     * @param[in]  algo_id          The protocol identifier override
     */
    void write_d2o_file(const std::string& filename,
                        const std::vector<fpt>& gridptr,
                        std::array<unsigned int, 3>& grid_dimensions,
                        const MatrixUnitcell& mat,
                        D2OFormat::CompressionAlgo algo_id = D2OFormat::CompressionAlgo::AUTO);

    /**
     * @brief      Compress a stream
     *
     * @param[in]  originstr  Original string to compress
     * @param[in]  algo_id    The algorithm identifier
     *
     * @return     Compressed stream
     */
    std::string compress_stream(const std::string& originstr, D2OFormat::CompressionAlgo algo_id);

    /**
     * @brief      Compress data stream using all possible compression algos
     *
     * @param[in]  origin  Uncompressed data stream
     *
     * @return     Unordered map containing compressed streams as strings
     */
    std::unordered_map<D2OFormat::CompressionAlgo, std::string, std::hash<D2OFormat::CompressionAlgo>> d2o_compress_all(const std::string& originstr);

    /**
     * @brief      Check that a compressed stream can be correctly decompressed
     *
     * @param[in]  algo_id          Which compression algorithm has been used
     * @param[in]  verificationstr  Original data stream for verification (string)
     * @param[in]  compressedstr    Compressed data stream (string)
     *
     * @return     Whether compressed stream can be correctly decompressed
     */
    bool check_decompression(D2OFormat::CompressionAlgo algo_id, const std::string& verificationstr, const std::string& compressedstr);

}

#endif // _D2O_FORMAT