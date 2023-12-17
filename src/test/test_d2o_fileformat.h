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

#ifndef _TEST_D2O_FileFormat
#define _TEST_D2O_FileFormat

#include <cppunit/extensions/HelperMacros.h>
#include <random>

#include "scalar_field.h"

class TestD2OFileFormat : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TestD2OFileFormat );
  CPPUNIT_TEST( test_gzip_compression );
  CPPUNIT_TEST( test_lzma_compression );
  CPPUNIT_TEST( test_bzip2_compression );
  CPPUNIT_TEST( test_autocompression );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void test_gzip_compression();
  void test_lzma_compression();
  void test_bzip2_compression();
  void test_autocompression();

private:
  uint32_t get_protocol_id(const std::string& filename);
};

#endif  // _TEST_D2O_FileFormat