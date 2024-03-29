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

#include "scalar_field.h"
#include "generator.h"

class TestGenerator : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TestGenerator );
  CPPUNIT_TEST( test_generator_genus );
  CPPUNIT_TEST( test_generator_benzene );
  CPPUNIT_TEST_SUITE_END();

public:
    void setUp();
    void test_generator_genus();
    void test_generator_benzene();
    void tearDown();

private:
};

#endif  // _TEST_ISOSURFACE
