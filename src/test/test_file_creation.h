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

#ifndef _TEST_FILE_CREATION
#define _TEST_FILE_CREATION

#include <cppunit/extensions/HelperMacros.h>

// we need these for the MD5 checksums
#include <iomanip>
#include <openssl/md5.h>
#include <openssl/evp.h>

#include "scalar_field.h"
#include "isosurface.h"
#include "isosurface_mesh.h"

class TestFileCreation : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE( TestFileCreation );
    CPPUNIT_TEST( test_ply_file );
    CPPUNIT_TEST( test_stl_file );
    CPPUNIT_TEST( test_obj_file );
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp();
    void test_ply_file();
    void test_stl_file();
    void test_obj_file();
    void tearDown();

private:
    std::string md5(const std::string &str);
};

#endif