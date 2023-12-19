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

#ifndef _TEST_ISOSURFACE
#define _TEST_ISOSURFACE

#include <cppunit/extensions/HelperMacros.h>
#include <memory>
#include <unistd.h>
#include <filesystem>

#include "scalar_field.h"
#include "isosurface.h"
#include "isosurface_mesh.h"

class TestIsosurface : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TestIsosurface );
  CPPUNIT_TEST( test_marching_cubes );
  CPPUNIT_TEST( test_marching_tetrahedra );
  CPPUNIT_TEST( test_obj );
  CPPUNIT_TEST( test_ply );
  CPPUNIT_TEST( test_stl );
  CPPUNIT_TEST_SUITE_END();

public:
    void setUp();
    void test_marching_cubes();
    void test_marching_tetrahedra();
    void test_obj();
    void test_ply();
    void test_stl();
    void tearDown();

private:
    std::unique_ptr<ScalarField> sf;
};

#endif  // _TEST_ISOSURFACE
