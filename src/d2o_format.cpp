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
void D2OFormat::write_d2o_file(const std::string& filename,
                               const std::vector<fpt>& gridptr,
                               std::array<unsigned int, 3>& grid_dimensions,
                               const MatrixUnitcell& mat,
                               uint32_t protocol_override) {
    // write file format token
    std::ofstream outfile(filename, std::ios::binary);
    char buf[] = "D2O";
    outfile.write(buf, 3);

    // auto-look for best compression algo
    std::cout << "Looking for best compression algorithm." << std::endl;

    // copy data to char array
    size_t gridptrsz = gridptr.size() * sizeof(fpt);
    char* data = new char[gridptrsz];
    memcpy(data, &gridptr[0], gridptrsz);
    const std::string originstr(data, gridptrsz);

    // determine best compression format
    std::vector<std::string> compstr = d2o_compress_all(originstr);
    std::vector<unsigned int> stringsize(3,0);

    // clean up data object (no longer needed)
    delete[] data;

    // list of compression algos
    static const std::vector<std::string> algos = {
        "GZIP",
        "LZMA",
        "BZIP2"
    };


    unsigned int idx = 0;
    if(protocol_override == 0) {
        // determine size for each type of compression
        for(unsigned int i=0; i<3; i++) {
            stringsize[i] = compstr[i].size();
        }

        // find best compression algo and base protocol id on this
        idx = std::distance(std::begin(stringsize), std::min_element(std::begin(stringsize), std::end(stringsize)));
        std::cout << "Best algorithm: " << algos[idx] << std::endl;
    } else {
        idx = protocol_override - 1;
        std::cout << "Overruling compression algo to: " << algos[idx] << std::endl;
    }

    // verify that the compressed stream can be correctly decompressed
    bool valid_decompression = check_decompression(idx, originstr, compstr[idx]);
    if(!valid_decompression) {
        throw std::runtime_error("Decompression could not be verified. Please try again or force to use a different compression algorithm.");
    } else {
        std::cout << "Decompression verification passed." << std::endl;
    }

    // capture compressed data
    uint32_t protocol_id = idx + 1;

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

    // write compressed data size and compressed data
    uint64_t sz = compstr[idx].size();
    outfile.write((char*)&sz, sizeof(uint64_t));
    outfile.write(compstr[idx].data(), compstr[idx].size());

    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    std::cout << "Writing " << filename << " ("
    << (boost::format("%0.1f") % ((float)size / 1024.f)).str()
    << "kb)." << std::endl;
}

std::vector<std::string> D2OFormat::d2o_compress_all(const std::string& originstr) {
    std::vector<std::string> compressed_strings;
    unsigned int gridptrsz = originstr.size();

    for(unsigned int i=0; i<3; i++) {
        std::istringstream origin(originstr);
        boost::iostreams::filtering_istreambuf in;
        switch(i) {
            case 0:
            {
                std::cout << "Trying GZIP: ";
                in.push(
                    boost::iostreams::gzip_compressor(
                        boost::iostreams::gzip_params(
                            boost::iostreams::gzip::best_compression
                        )
                    )
                );
            }
            break;
            case 1:
            {
                std::cout << "Trying LZMA: ";
                in.push(
                    boost::iostreams::lzma_compressor(
                        boost::iostreams::lzma_params(
                            boost::iostreams::lzma::best_compression
                        )
                    )
                );
            }
            break;
            case 2:
            {
                std::cout << "Trying BZIP2: ";
                in.push(
                    boost::iostreams::bzip2_compressor()
                );
            }
            break;
        }

        in.push(origin);
        std::ostringstream compressed;
        boost::iostreams::copy(in, compressed);
        compressed_strings.emplace_back(compressed.str());
        size_t sz = compressed_strings.back().size();

        std::cout << (boost::format("%0.1f") % ((float)sz / 1024.f)).str()
        << " kb ("
        << (boost::format("%0.2f") % ((float)sz / (float)gridptrsz * 100.0f)).str()
        << " %)." << std::endl;
    }

    return compressed_strings;
}

/**
 * @brief      Check that a compressed stream can be correctly decompressed
 *
 * @param[in]  protocol         Which compression algorithm has been used
 * @param[in]  verificationstr  Original data stream (string)
 * @param[in]  compressedstr    Compressed data stream (string)
 *
 * @return     Whether compressed stream can be correctly decompressed
 */
bool D2OFormat::check_decompression(uint32_t protocol, const std::string& verificationstr, const std::string& compressedstr) {
    // decompress
    std::istringstream compressed(compressedstr);
    boost::iostreams::filtering_istreambuf in;

    switch(protocol) {
        case 0:
        {
            // GZIP compression
            std::cout << "Building GZIP decompressor" << std::endl;
            in.push(boost::iostreams::gzip_decompressor());
        }
        break;
        case 1:
        {
            std::cout << "Building LZMA decompressor" << std::endl;
            in.push(boost::iostreams::lzma_decompressor());
        }
        break;
        case 2:
        {
            std::cout << "Building BZIP2 decompressor" << std::endl;
            in.push(boost::iostreams::bzip2_decompressor());
        }
        break;
        default:
            throw std::runtime_error("Invalid algorithm id encountered for decompression verification: " + std::to_string(protocol));
        break;
    }

    in.push(compressed);

    std::ostringstream origin;
    boost::iostreams::copy(in, origin);
    const std::string originstr = origin.str();

    return verificationstr == originstr;
}