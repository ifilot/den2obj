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
                               D2OFormat::CompressionAlgo algo_id) {
    // write file format token
    std::ofstream outfile(filename, std::ios::binary);
    char buf[] = "D2O";
    outfile.write(buf, 3);

    std::cout << "Writing .D2O file, using compression algo id: " << (int)algo_id << std::endl;

    // copy data to char array
    size_t gridptrsz = gridptr.size() * sizeof(fpt);
    char* data = new char[gridptrsz];
    memcpy(data, &gridptr[0], gridptrsz);
    const std::string originstr(data, gridptrsz);

    std::string compressedstr;
    if(algo_id == D2OFormat::CompressionAlgo::AUTO) { // check which compression algo works best
        // auto-look for best compression algo
        std::cout << "Looking for best compression algorithm." << std::endl;

        // determine best compression format
        auto compstr = d2o_compress_all(originstr);
        std::unordered_map<D2OFormat::CompressionAlgo, unsigned int, std::hash<D2OFormat::CompressionAlgo>> stringsizes;

        // determine size for each type of compression
        for(const auto& i : compstr) {
            stringsizes.emplace(i.first, i.second.size());
        }

        // find best compression algo and base protocol id on this
        auto got = std::min_element(stringsizes.begin(), stringsizes.end(), [](const auto& l, const auto& r) { return l.second < r.second; });
        compressedstr = compstr.find(got->first)->second;

        // overwrite algo_id
        algo_id = got->first;
        // std::cout << "Best algorithm: " << D2OFormat::algos[idx] << std::endl;
    } else { // use the user-supplied compression algo
        compressedstr = D2OFormat::compress_stream(originstr, algo_id);
        // std::cout << "Overruling compression algo to: " << D2OFormat::algos[idx] << std::endl;
    }

    // verify that the compressed stream can be correctly decompressed
    bool valid_decompression = check_decompression(algo_id, originstr, compressedstr);
    if(!valid_decompression) {
        throw std::runtime_error("Decompression could not be verified. Please try again or force to use a different compression algorithm.");
    } else {
        std::cout << "Decompression verification passed." << std::endl;
    }

    // capture compressed data
    uint32_t protocol_id = static_cast<uint32_t>(algo_id);

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
    uint64_t sz = compressedstr.size();
    outfile.write((char*)&sz, sizeof(uint64_t));
    outfile.write(compressedstr.data(), compressedstr.size());

    // clean up data object (no longer needed)
    delete[] data;

    // close file
    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    std::cout << "Writing " << filename << " ("
    << (boost::format("%0.1f") % ((float)size / 1024.f)).str()
    << "kb)." << std::endl;
}

std::unordered_map<D2OFormat::CompressionAlgo, std::string, std::hash<D2OFormat::CompressionAlgo>> D2OFormat::d2o_compress_all(const std::string& originstr) {
    std::unordered_map<D2OFormat::CompressionAlgo, std::string, std::hash<D2OFormat::CompressionAlgo>> compressed_strings;

    for(const auto& i : D2OFormat::algos) {
        // skip auto tag
        if(i.second == D2OFormat::CompressionAlgo::AUTO) {
            continue;
        }

        std::cout << "Trying " << i.first << ": ";
        
        const std::string compressedstr = D2OFormat::compress_stream(originstr, i.second);
        compressed_strings.emplace(i.second, compressedstr);
        size_t sz = compressedstr.size();

        std::cout << (boost::format("%0.1f") % ((float)sz / 1024.f)).str()
        << " kb ("
        << (boost::format("%0.2f") % ((float)sz / (float)originstr.size() * 100.0f)).str()
        << " %)." << std::endl;

        if(sz == 0) {
            throw std::runtime_error("A compression size of 0 is incorrect. Something is wrong with this file.");
        }
    }

    return compressed_strings;
}

/**
 * @brief      Compress a stream
 *
 * @param[in]  originstr  Original string to compress
 * @param[in]  algo_id    The algorithm identifier
 *
 * @return     Compressed stream
 */
std::string D2OFormat::compress_stream(const std::string& originstr, D2OFormat::CompressionAlgo algo_id) {
    std::istringstream origin(originstr);
    boost::iostreams::filtering_istreambuf in;

    switch(algo_id) {
        case D2OFormat::CompressionAlgo::GZIP:
        {
            in.push(
                boost::iostreams::gzip_compressor(
                    boost::iostreams::gzip_params(
                        boost::iostreams::gzip::best_compression
                    )
                )
            );
        }
        break;
        case D2OFormat::CompressionAlgo::LZMA:
        {
            in.push(
                boost::iostreams::lzma_compressor(
                    boost::iostreams::lzma_params(
                        boost::iostreams::lzma::default_compression
                    )
                )
            );
        }
        break;
        case D2OFormat::CompressionAlgo::BZIP2:
        {
            in.push(
                boost::iostreams::bzip2_compressor()
            );
        }
        break;
        default:
            throw std::runtime_error("Invalid algorithm id received.");
        break;
    }

    in.push(origin);
    std::ostringstream compressed;
    boost::iostreams::copy(in, compressed);
    return compressed.str();
}

/**
 * @brief      Check that a compressed stream can be correctly decompressed
 *
 * @param[in]  algo_id          Which compression algorithm has been used
 * @param[in]  verificationstr  Original data stream (string)
 * @param[in]  compressedstr    Compressed data stream (string)
 *
 * @return     Whether compressed stream can be correctly decompressed
 */
bool D2OFormat::check_decompression(D2OFormat::CompressionAlgo algo_id, const std::string& verificationstr, const std::string& compressedstr) {
    // decompress
    std::istringstream compressed(compressedstr);
    boost::iostreams::filtering_istreambuf in;

    switch(algo_id) {
        case D2OFormat::CompressionAlgo::GZIP:
        {
            // GZIP compression
            std::cout << "Building GZIP decompressor" << std::endl;
            in.push(boost::iostreams::gzip_decompressor());
        }
        break;
        case D2OFormat::CompressionAlgo::LZMA:
        {
            std::cout << "Building LZMA decompressor" << std::endl;
            in.push(boost::iostreams::lzma_decompressor());
        }
        break;
        case D2OFormat::CompressionAlgo::BZIP2:
        {
            std::cout << "Building BZIP2 decompressor" << std::endl;
            in.push(boost::iostreams::bzip2_decompressor());
        }
        break;
        default:
            throw std::runtime_error("Invalid algorithm id encountered for decompression verification.");
        break;
    }

    in.push(compressed);

    std::ostringstream origin;
    boost::iostreams::copy(in, origin);
    const std::string originstr = origin.str();

    return verificationstr == originstr;
}