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

#ifndef _TEST_SCALARFIELD
#define _TEST_SCALARFIELD

#include <cppunit/extensions/HelperMacros.h>

#include "scalar_field.h"

class TestScalarField : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TestScalarField );
  CPPUNIT_TEST( testReadingCHGCAR );
  CPPUNIT_TEST( testReadingCUB );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void testReadingCHGCAR();
  void testReadingCUB();

private:
};

#endif  // _TEST_SCALARFIELD
