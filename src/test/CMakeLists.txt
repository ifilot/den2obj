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

if(STATIC_DEPENDENCIES)
    find_path(CPPUNIT_INCLUDE_DIR NAMES cppunit/TestFixture.h REQUIRED)
    find_library(CPPUNIT_LIBRARY NAMES cppunit)
    if(NOT CPPUNIT_LIBRARY)
        message(FATAL_ERROR "Could not locate static dependency CPPUNIT_LIBRARY")
    endif()
    add_library(Den2ObjCppUnit UNKNOWN IMPORTED)
    set_target_properties(Den2ObjCppUnit PROPERTIES
        IMPORTED_LOCATION "${CPPUNIT_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CPPUNIT_INCLUDE_DIR}")
    set(CPPUNIT_TARGET Den2ObjCppUnit ${CMAKE_DL_LIBS})
else()
    pkg_check_modules(CPPUNIT REQUIRED IMPORTED_TARGET cppunit)
    set(CPPUNIT_TARGET PkgConfig::CPPUNIT)
endif()

add_library(unittest STATIC unittest.cpp)
target_compile_features(unittest PUBLIC cxx_std_17)
target_link_libraries(unittest PRIVATE ${CPPUNIT_TARGET})

set(EXECUTABLES TestIsosurface TestScalarField TestD2OFileFormat TestGenerator)

if(USE_GCOV)
    list(APPEND EXECUTABLES TestFileCreation)
endif()

add_executable(TestIsosurface test_isosurface.cpp)
add_executable(TestScalarField test_scalarfield.cpp)
add_executable(TestD2OFileFormat test_d2o_fileformat.cpp)
add_executable(TestGenerator test_generator.cpp)

if(USE_GCOV)
    add_executable(TestFileCreation test_file_creation.cpp)
endif()

add_test(NAME DatasetSetup COMMAND tar --no-same-owner -xvjf dataset.tar.bz2)
add_test(NAME DatasetCleanup COMMAND rm -rvf CHGCAR_* PARCHG_* *.d2o *.cub)
set_tests_properties(DatasetSetup PROPERTIES FIXTURES_SETUP Dataset)
set_tests_properties(DatasetCleanup PROPERTIES FIXTURES_CLEANUP Dataset)

foreach(testexec ${EXECUTABLES})
    target_compile_features(${testexec} PUBLIC cxx_std_17)
    if(USE_GCOV)
        target_compile_options(${testexec} PRIVATE --coverage)
        target_link_options(${testexec} PRIVATE --coverage)
    endif()

    target_link_libraries(${testexec}
        PRIVATE
            unittest
            den2objsources
            ${CPPUNIT_TARGET}
    )

    add_test(NAME ${testexec} COMMAND ${testexec})
    set_tests_properties(${testexec} PROPERTIES FIXTURES_REQUIRED Dataset)
endforeach()

configure_file(testinput/dataset.tar.bz2 dataset.tar.bz2 COPYONLY)
