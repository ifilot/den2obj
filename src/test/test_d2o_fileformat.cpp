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

void TestD2OFileFormat::setUp() {
    ScalarField sf("CHGCAR_CH4", ScalarFieldInputFileType::SFF_CHGCAR);
    CPPUNIT_ASSERT_EQUAL( (uint)0, sf.get_size() );

    static const std::string filename = "chgcar_ch4_base.d2o";
    sf.read();
    sf.write_d2o_binary(filename);
}

void TestD2OFileFormat::test_gzip_compression() {
    // create scalar field
    ScalarField sf("chgcar_ch4_base.d2o", ScalarFieldInputFileType::SFF_D2O);

    static const std::string filename = "chgcar_ch4_gzip.d2o";

    sf.write_d2o_binary(filename, 1);
    CPPUNIT_ASSERT_EQUAL((uint32_t)1, this->get_protocol_id(filename));

    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare vector points
    for(unsigned int i=0; i<grid.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[i], sf.get_grid()[i], 1e-12);
    }
}

void TestD2OFileFormat::test_lzma_compression() {
    // create scalar field
    ScalarField sf("chgcar_ch4_base.d2o", ScalarFieldInputFileType::SFF_CHGCAR);

    static const std::string filename = "chgcar_ch4_lzma.d2o";

    sf.write_d2o_binary(filename, 2);
    CPPUNIT_ASSERT_EQUAL((uint32_t)2, this->get_protocol_id(filename));

    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare vector points
    for(unsigned int i=0; i<grid.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[i], sf.get_grid()[i], 1e-12);
    }
}

void TestD2OFileFormat::test_bzip2_compression() {
    // create scalar field
    ScalarField sf("chgcar_ch4_base.d2o", ScalarFieldInputFileType::SFF_CHGCAR);

    static const std::string filename = "chgcar_ch4_bzip2.d2o";

    sf.write_d2o_binary(filename, 3);
    CPPUNIT_ASSERT_EQUAL((uint32_t)3, this->get_protocol_id(filename));

    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare vector points
    for(unsigned int i=0; i<grid.size(); i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[i], sf.get_grid()[i], 1e-12);
    }
}

void TestD2OFileFormat::test_autocompression() {
    // create scalar field
    ScalarField sf("chgcar_ch4_base.d2o", ScalarFieldInputFileType::SFF_CHGCAR);
    CPPUNIT_ASSERT_EQUAL( (uint)0, sf.get_size() );

    static const std::string filename = "chgcar_ch4_auto.d2o";

    sf.write_d2o_binary(filename, 0);
    CPPUNIT_ASSERT_EQUAL((uint32_t)2, this->get_protocol_id(filename));

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