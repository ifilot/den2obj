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

#include "test_generator.h"

CPPUNIT_TEST_SUITE_REGISTRATION( TestGenerator );

void TestGenerator::setUp() {
}

void TestGenerator::tearDown() {
}

void TestGenerator::test_generator_genus() {
    Generator gen;
    const std::string filename = "genus2.d2o";

    // generator dataset
    gen.build_dataset("genus2", filename, D2OFormat::CompressionAlgo::BZIP2);

    // check that it can be read
    CPPUNIT_ASSERT_NO_THROW(ScalarField(filename, ScalarFieldInputFileType::SFF_D2O));
}

void TestGenerator::test_generator_benzene() {
    Generator gen;
    const std::string filename = "benzene_homo.d2o";

    // generator dataset
    gen.build_dataset("benzene_homo", filename, D2OFormat::CompressionAlgo::BZIP2);

    // check that it can be read
    CPPUNIT_ASSERT_NO_THROW(ScalarField(filename, ScalarFieldInputFileType::SFF_D2O));
}
