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

#include "d2o_format.h"

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
                    const MatrixUnitcell& mat) {
    // write file format token
    std::ofstream outfile(filename, std::ios::binary);
    char buf[] = "D2O";
    outfile.write(buf, 3);

    if(protocol_id == 0 || protocol_id > 2) {
        throw std::runtime_error("Invalid protocol id for d2o file: " + std::to_string(protocol_id));
    }

    // write protocol id
    outfile.write((char*)&protocol_id, sizeof(uint32_t));

    // write unit cell matrix
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            float val = mat(i,j);
            outfile.write((char*)&val, sizeof(float));
        }
    }

    // write number of data points
    for(unsigned int i=0; i<3; i++) {
        uint32_t sz = grid_dimensions[i];
        outfile.write((char*)&sz, sizeof(uint32_t));
    }

    // write floating point size
    uint8_t fptsz = sizeof(fpt);
    outfile.write((char*)&fptsz, sizeof(uint8_t));

    std::cout << "Floating point size determined at: " << (int)fptsz << " bytes" << std::endl;

    // compressing stream
    size_t gridptrsz = gridptr.size() * sizeof(fpt);
    char* data = new char[gridptrsz];
    memcpy(data, &gridptr[0], gridptrsz);
    std::istringstream origin(std::string(data, gridptrsz));
    boost::iostreams::filtering_istreambuf in;

    switch(protocol_id) {
        case 1:
            // GZIP compression
            in.push(
                boost::iostreams::gzip_compressor(
                    boost::iostreams::gzip_params(
                        boost::iostreams::gzip::best_compression
                    )
                )
            );
            std::cout << "Using GZIP compression." << std::endl;
        break;
        case 2:
            in.push(
                boost::iostreams::lzma_compressor()
            );
            std::cout << "Using LZMA compression." << std::endl;
        break;
        default:
            throw std::runtime_error("Invalid protocol id for d2o file: " + std::to_string(protocol_id));
        break;
    }

    in.push(origin);

    // store compression
    std::ostringstream compressed;
    boost::iostreams::copy(in, compressed);

    // output to file
    const std::string griddata = compressed.str();
    uint64_t sz = griddata.size();
    std::cout << "Compressed data to "
        << (boost::format("%0.1f") % ((float)sz / 1024.f)).str()
        << " kb ("
        << (boost::format("%0.2f") % ((float)sz / (float)gridptrsz * 100.0f)).str()
        << " %)." << std::endl;

    // write data size and data
    outfile.write((char*)&sz, sizeof(uint64_t));
    outfile.write(griddata.data(), griddata.size());

    // clean up
    delete[] data;
    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    std::cout << "Writing " << filename << " ("
    << (boost::format("%0.1f") % ((float)size / 1024.f)).str()
    << "kb)." << std::endl;
}
