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

#include "test_d2o_fileformat.h"

CPPUNIT_TEST_SUITE_REGISTRATION( TestD2OFileFormat );

void TestD2OFileFormat::setUp() {
    ScalarField sf("CHGCAR_CH4", ScalarFieldInputFileType::SFF_CHGCAR);
    CPPUNIT_ASSERT_EQUAL( (uint)0, sf.get_size() );

    static const std::string filename = basefile;
    sf.read();
    sf.write_d2o_binary(filename);
}

void TestD2OFileFormat::test_gzip_compression() {
    this->assert_compression_roundtrip(D2OFormat::CompressionAlgo::GZIP, 1, "chgcar_ch4_gzip.d2o");
}

void TestD2OFileFormat::test_lzma_compression() {
    this->assert_compression_roundtrip(D2OFormat::CompressionAlgo::LZMA, 2, "chgcar_ch4_lzma.d2o");
}

void TestD2OFileFormat::test_bzip2_compression() {
    this->assert_compression_roundtrip(D2OFormat::CompressionAlgo::BZIP2, 3, "chgcar_ch4_bzip2.d2o");
}

void TestD2OFileFormat::test_zstd_compression() {
    this->assert_compression_roundtrip(D2OFormat::CompressionAlgo::ZSTD, 4, "chgcar_ch4_zstd.d2o");
}

void TestD2OFileFormat::test_blosc_compression() {
    this->assert_compression_roundtrip(D2OFormat::CompressionAlgo::BLOSC, 5, "chgcar_ch4_blosc.d2o");
}

void TestD2OFileFormat::test_autocompression() {
    // create scalar field
    ScalarField sf(basefile, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL( (uint)1000000, sf.get_size() );

    static const std::string filename = "chgcar_ch4_auto.d2o";

    sf.write_d2o_binary(filename, D2OFormat::CompressionAlgo::AUTO);
    CPPUNIT_ASSERT_GREATEREQUAL((uint32_t)1, this->get_protocol_id(filename));
    CPPUNIT_ASSERT_LESSEQUAL((uint32_t)5, this->get_protocol_id(filename));

    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare vector points
    for(unsigned int i=0; i<grid.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[i], sf.get_grid()[i], 1e-12);
    }
}

void TestD2OFileFormat::assert_compression_roundtrip(D2OFormat::CompressionAlgo algo_id,
                                                     uint32_t protocol_id,
                                                     const std::string& filename) {
    // create scalar field
    ScalarField sf(basefile, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL( (uint)1000000, sf.get_size() );

    sf.write_d2o_binary(filename, algo_id);
    CPPUNIT_ASSERT_EQUAL(protocol_id, this->get_protocol_id(filename));

    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare vector points
    for(unsigned int i=0; i<grid.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[i], sf.get_grid()[i], 1e-12);
    }
}

uint32_t TestD2OFileFormat::get_protocol_id(const std::string& filename) {
    char buf[3];
    std::ifstream f(filename);
    f.read(buf, 3);
    uint32_t protocol_id = 0;
    f.read((char*)&protocol_id, sizeof(uint32_t));
    f.close();

    return protocol_id;
}

void TestD2OFileFormat::tearDown() {
    // do nothing
}
