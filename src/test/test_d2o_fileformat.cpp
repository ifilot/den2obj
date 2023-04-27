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

void TestD2OFileFormat::test_gzip_compression() {
    // create scalar field
    ScalarField sf("CHGCAR_CH4", ScalarFieldInputFileType::SFF_CHGCAR);
    CPPUNIT_ASSERT_EQUAL( (uint)0, sf.get_size() );

    // read atoms and check this
    sf.read_header_and_atoms();
    sf.read();

    static const std::string filename = "chgcar_ch4_gzip.d2o";

    sf.write_d2o_binary(filename, 1);
    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare random points
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,nrgridpts-1);
    for(unsigned int i=0; i<1000; i++) {
        unsigned int idx = dist(rng);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[idx], sf.get_grid()[idx], 1e-12);
    }
}

void TestD2OFileFormat::test_lzma_compression() {
    // create scalar field
    ScalarField sf("CHGCAR_CH4", ScalarFieldInputFileType::SFF_CHGCAR);
    CPPUNIT_ASSERT_EQUAL( (uint)0, sf.get_size() );

    // read atoms and check this
    sf.read_header_and_atoms();
    sf.read();

    static const std::string filename = "chgcar_ch4_lzma.d2o";

    sf.write_d2o_binary(filename, 2);
    const auto nrgridpts = sf.get_grid().size();
    const auto grid = sf.get_grid();

    sf = ScalarField(filename, ScalarFieldInputFileType::SFF_D2O);
    CPPUNIT_ASSERT_EQUAL(nrgridpts, sf.get_grid().size());

    // compare random points
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,nrgridpts-1);
    for(unsigned int i=0; i<1000; i++) {
        unsigned int idx = dist(rng);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grid[idx], sf.get_grid()[idx], 1e-12);
    }
}
