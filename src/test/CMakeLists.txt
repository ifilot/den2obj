 #*************************************************************************
 #   CMakeLists.txt  --  This file is part of DEN2OBJ.                    *
 #                                                                        *
 #   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 #                                                                        *
 #   DEN2OBJ is free software: you can redistribute it and/or modify      *
 #   it under the terms of the GNU General Public License as published    *
 #   by the Free Software Foundation, either version 3 of the License,    *
 #   or (at your option) any later version.                               *
 #                                                                        *
 #   DEN2OBJ is distributed in the hope that it will be useful,           *
 #   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 #   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 #   See the GNU General Public License for more details.                 *
 #                                                                        *
 #   You should have received a copy of the GNU General Public License    *
 #   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 #                                                                        *
 #*************************************************************************/

add_library(unittest STATIC unittest.cpp)

# set executables
SET(EXECUTABLES TestIsosurface TestScalarField TestD2OFileFormat)

#######################################################
# Add executables
#######################################################
add_executable(TestIsosurface test_isosurface.cpp)
add_executable(TestScalarField test_scalarfield.cpp)
add_executable(TestD2OFileFormat test_d2o_fileformat.cpp)

#######################################################
# Link mkmsources and other dependencies
#######################################################
if(CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CPPUNIT_LIB ${CPPUNIT_DEBUG_LIBRARY})
else()
    SET(CPPUNIT_LIB ${CPPUNIT_LIBRARY})
endif()

# add unpacking of dataset to the test suite
add_test(NAME DatasetSetup COMMAND tar -xvjf dataset.tar.bz2)
add_test(NAME DatasetCleanup COMMAND rm -rvf CHGCAR_* PARCHG_*)
set_tests_properties(DatasetSetup PROPERTIES FIXTURES_SETUP Dataset)
set_tests_properties(DatasetCleanup PROPERTIES FIXTURES_CLEANUP Dataset)

# configure common settings for test executables
foreach(testexec ${EXECUTABLES})
    target_link_libraries(${testexec}
                          unittest
                          den2objsources
                          ${OPENVDB_LIBS}
                          ${Boost_LIBRARIES}
                          ${CPPUNIT_LIB}
                          ${LIBLZMA_LIBRARIES}
                          ${ZLIB_LIBRARIES}
                          -lgcov)
    set_target_properties(${testexec} PROPERTIES COMPILE_FLAGS "--coverage")
    add_test(NAME ${testexec} COMMAND ${testexec})
    set_tests_properties(${testexec} PROPERTIES FIXTURES_REQUIRED Dataset)
endforeach()

#######################################################
# add testinput files
#######################################################
configure_file(testinput/dataset.tar.bz2 dataset.tar.bz2 COPYONLY)
