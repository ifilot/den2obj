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
#include "cli_format.h"

#include <blosc.h>
#include <iomanip>
#include <zstd.h>

namespace {
    std::vector<std::pair<D2OFormat::CompressionAlgo, std::string>> ordered_compression_algos() {
        return {
            {D2OFormat::CompressionAlgo::GZIP, "gzip"},
            {D2OFormat::CompressionAlgo::LZMA, "lzma"},
            {D2OFormat::CompressionAlgo::BZIP2, "bzip2"},
            {D2OFormat::CompressionAlgo::ZSTD, "zstd"},
            {D2OFormat::CompressionAlgo::BLOSC, "blosc"}
        };
    }

    void print_compression_table(
        const std::unordered_map<D2OFormat::CompressionAlgo, std::string, std::hash<D2OFormat::CompressionAlgo>>& compressed_strings,
        D2OFormat::CompressionAlgo selected_algo,
        size_t original_size) {
        CLIFormat::print_section("Compression benchmark");
        std::cout << "+--------+-----------+--------+----------+" << std::endl;
        std::cout << "| Method | Size      | Ratio  | Status   |" << std::endl;
        std::cout << "+--------+-----------+--------+----------+" << std::endl;

        for(const auto& algo : ordered_compression_algos()) {
            const auto got = compressed_strings.find(algo.first);
            if(got == compressed_strings.end()) {
                continue;
            }

            const double ratio = static_cast<double>(got->second.size()) /
                                 static_cast<double>(original_size) * 100.0;

            const std::string status = algo.first == selected_algo ? "best" : "";
            const std::string padded_status = status + std::string(8 - status.size(), ' ');

            std::cout << "| " << std::left << std::setw(6) << algo.second
                      << " | " << std::right << std::setw(9) << CLIFormat::format_kb(got->second.size())
                      << " | " << std::setw(5) << (boost::format("%0.2f") % ratio).str() << "%"
                      << " | " << CLIFormat::colorize(padded_status, CLIFormat::Color::GREEN)
                      << " |" << std::endl;
        }

        std::cout << "+--------+-----------+--------+----------+" << std::endl;
        CLIFormat::print_kv("Selected", D2OFormat::compression_algo_name(selected_algo));
    }
}

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

    CLIFormat::print_section("Writing D2O file");
    CLIFormat::print_kv("Compression", D2OFormat::compression_algo_name(algo_id));

    // copy data to char array
    size_t gridptrsz = gridptr.size() * sizeof(fpt);
    char* data = new char[gridptrsz];
    memcpy(data, &gridptr[0], gridptrsz);
    const std::string originstr(data, gridptrsz);

    std::string compressedstr;
    if(algo_id == D2OFormat::CompressionAlgo::AUTO) { // check which compression algo works best
        // determine best compression format
        auto compstr = d2o_compress_all(originstr);
        std::unordered_map<D2OFormat::CompressionAlgo, size_t, std::hash<D2OFormat::CompressionAlgo>> stringsizes;

        // determine size for each type of compression
        for(const auto& i : compstr) {
            stringsizes.emplace(i.first, i.second.size());
        }

        // find best compression algo and base protocol id on this
        auto got = std::min_element(stringsizes.begin(), stringsizes.end(), [](const auto& l, const auto& r) { return l.second < r.second; });
        compressedstr = compstr.find(got->first)->second;

        // overwrite algo_id
        algo_id = got->first;
        print_compression_table(compstr, algo_id, originstr.size());
    } else { // use the user-supplied compression algo
        compressedstr = D2OFormat::compress_stream(originstr, algo_id);
    }

    // verify that the compressed stream can be correctly decompressed
    bool valid_decompression = check_decompression(algo_id, originstr, compressedstr);
    if(!valid_decompression) {
        throw std::runtime_error("Decompression could not be verified. Please try again or force to use a different compression algorithm.");
    } else {
        CLIFormat::print_kv("Verification", "passed");
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
    CLIFormat::print_kv("Float size", std::to_string((int)fptsz) + " bytes");

    // write compressed data size and compressed data
    uint64_t sz = compressedstr.size();
    outfile.write((char*)&sz, sizeof(uint64_t));
    outfile.write(compressedstr.data(), compressedstr.size());

    // clean up data object (no longer needed)
    delete[] data;

    // close file
    outfile.close();

    std::uintmax_t size = boost::filesystem::file_size(filename);
    CLIFormat::print_kv("Output", filename);
    CLIFormat::print_kv("File size", CLIFormat::format_kb(size));
}

std::unordered_map<D2OFormat::CompressionAlgo, std::string, std::hash<D2OFormat::CompressionAlgo>> D2OFormat::d2o_compress_all(const std::string& originstr) {
    std::unordered_map<D2OFormat::CompressionAlgo, std::string, std::hash<D2OFormat::CompressionAlgo>> compressed_strings;

    for(const auto& i : ordered_compression_algos()) {
        const std::string compressedstr = D2OFormat::compress_stream(originstr, i.first);
        compressed_strings.emplace(i.first, compressedstr);
        size_t sz = compressedstr.size();

        if(sz == 0) {
            throw std::runtime_error("A compression size of 0 is incorrect. Something is wrong with this file.");
        }
    }

    return compressed_strings;
}

/**
 * @brief      Get a human-readable compression algorithm name
 *
 * @param[in]  algo_id  The algorithm identifier
 *
 * @return     Algorithm name
 */
std::string D2OFormat::compression_algo_name(D2OFormat::CompressionAlgo algo_id) {
    switch(algo_id) {
        case D2OFormat::CompressionAlgo::AUTO:
            return "auto";
        case D2OFormat::CompressionAlgo::GZIP:
            return "gzip";
        case D2OFormat::CompressionAlgo::LZMA:
            return "lzma";
        case D2OFormat::CompressionAlgo::BZIP2:
            return "bzip2";
        case D2OFormat::CompressionAlgo::ZSTD:
            return "zstd";
        case D2OFormat::CompressionAlgo::BLOSC:
            return "blosc";
        default:
            return "unknown";
    }
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
        case D2OFormat::CompressionAlgo::ZSTD:
        {
            const size_t compressed_bound = ZSTD_compressBound(originstr.size());
            if(ZSTD_isError(compressed_bound)) {
                throw std::runtime_error("Could not determine ZSTD compression bound.");
            }

            std::string compressed(compressed_bound, '\0');
            const size_t compressed_size = ZSTD_compress(&compressed[0],
                                                         compressed.size(),
                                                         originstr.data(),
                                                         originstr.size(),
                                                         19);
            if(ZSTD_isError(compressed_size)) {
                throw std::runtime_error("ZSTD compression failed: " + std::string(ZSTD_getErrorName(compressed_size)));
            }

            compressed.resize(compressed_size);
            return compressed;
        }
        break;
        case D2OFormat::CompressionAlgo::BLOSC:
        {
            if(originstr.size() > BLOSC_MAX_BUFFERSIZE) {
                throw std::runtime_error("Input is too large for Blosc compression.");
            }

            std::string compressed(originstr.size() + BLOSC_MAX_OVERHEAD, '\0');
            const int compressed_size = blosc_compress_ctx(9,
                                                           BLOSC_BITSHUFFLE,
                                                           sizeof(fpt),
                                                           originstr.size(),
                                                           originstr.data(),
                                                           &compressed[0],
                                                           compressed.size(),
                                                           "blosclz",
                                                           0,
                                                           1);
            if(compressed_size <= 0) {
                throw std::runtime_error("Blosc compression failed.");
            }

            compressed.resize(compressed_size);
            return compressed;
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
 * @brief      Decompress a stream
 *
 * @param[in]  algo_id        The algorithm identifier
 * @param[in]  compressedstr  Compressed data stream
 * @param[in]  expected_size  Expected uncompressed byte size; 0 disables validation
 *
 * @return     Decompressed stream
 */
std::string D2OFormat::decompress_stream(D2OFormat::CompressionAlgo algo_id,
                                         const std::string& compressedstr,
                                         size_t expected_size) {
    switch(algo_id) {
        case D2OFormat::CompressionAlgo::GZIP:
        case D2OFormat::CompressionAlgo::LZMA:
        case D2OFormat::CompressionAlgo::BZIP2:
        {
            std::istringstream compressed(compressedstr);
            boost::iostreams::filtering_istreambuf in;

            switch(algo_id) {
                case D2OFormat::CompressionAlgo::GZIP:
                    CLIFormat::print_kv("Decompressor", "gzip");
                    in.push(boost::iostreams::gzip_decompressor());
                break;
                case D2OFormat::CompressionAlgo::LZMA:
                    CLIFormat::print_kv("Decompressor", "lzma");
                    in.push(boost::iostreams::lzma_decompressor());
                break;
                case D2OFormat::CompressionAlgo::BZIP2:
                    CLIFormat::print_kv("Decompressor", "bzip2");
                    in.push(boost::iostreams::bzip2_decompressor());
                break;
                default:
                    break;
            }

            in.push(compressed);

            std::ostringstream origin;
            boost::iostreams::copy(in, origin);
            const std::string originstr = origin.str();

            if(expected_size != 0 && originstr.size() != expected_size) {
                throw std::runtime_error("Decompressed stream has an unexpected size.");
            }

            return originstr;
        }
        break;
        case D2OFormat::CompressionAlgo::ZSTD:
        {
            CLIFormat::print_kv("Decompressor", "zstd");

            const unsigned long long frame_size = ZSTD_getFrameContentSize(compressedstr.data(), compressedstr.size());
            if(frame_size == ZSTD_CONTENTSIZE_ERROR) {
                throw std::runtime_error("Invalid ZSTD compressed stream.");
            }

            size_t decompressed_size = expected_size;
            if(decompressed_size == 0) {
                if(frame_size == ZSTD_CONTENTSIZE_UNKNOWN) {
                    throw std::runtime_error("Cannot determine ZSTD decompressed size.");
                }
                decompressed_size = static_cast<size_t>(frame_size);
            } else if(frame_size != ZSTD_CONTENTSIZE_UNKNOWN && frame_size != decompressed_size) {
                throw std::runtime_error("ZSTD decompressed stream has an unexpected size.");
            }

            std::string originstr(decompressed_size, '\0');
            const size_t actual_size = ZSTD_decompress(&originstr[0],
                                                       originstr.size(),
                                                       compressedstr.data(),
                                                       compressedstr.size());
            if(ZSTD_isError(actual_size)) {
                throw std::runtime_error("ZSTD decompression failed: " + std::string(ZSTD_getErrorName(actual_size)));
            }

            if(actual_size != decompressed_size) {
                throw std::runtime_error("ZSTD decompressed stream has an unexpected size.");
            }

            return originstr;
        }
        break;
        case D2OFormat::CompressionAlgo::BLOSC:
        {
            CLIFormat::print_kv("Decompressor", "blosc");

            size_t decompressed_size = 0;
            size_t compressed_size = 0;
            size_t block_size = 0;
            blosc_cbuffer_sizes(compressedstr.data(), &decompressed_size, &compressed_size, &block_size);

            if(decompressed_size == 0 || compressed_size == 0) {
                throw std::runtime_error("Invalid Blosc compressed stream.");
            }

            if(expected_size != 0 && decompressed_size != expected_size) {
                throw std::runtime_error("Blosc decompressed stream has an unexpected size.");
            }

            std::string originstr(decompressed_size, '\0');
            const int actual_size = blosc_decompress_ctx(compressedstr.data(),
                                                         &originstr[0],
                                                         originstr.size(),
                                                         1);
            if(actual_size <= 0) {
                throw std::runtime_error("Blosc decompression failed.");
            }

            if(static_cast<size_t>(actual_size) != decompressed_size) {
                throw std::runtime_error("Blosc decompressed stream has an unexpected size.");
            }

            return originstr;
        }
        break;
        default:
            throw std::runtime_error("Invalid algorithm id encountered for decompression.");
        break;
    }
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
    const std::string originstr = D2OFormat::decompress_stream(algo_id, compressedstr, verificationstr.size());
    return verificationstr == originstr;
}
