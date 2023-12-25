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

# set minimum cmake requirements
cmake_minimum_required(VERSION 3.5)
project (den2obj)

# change compiler directives when code coverage is required
if(USE_GCOV)
  set(CMAKE_BUILD_TYPE "Debug")

  # Set global c and c++ flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")

  # Link flags used for creating executables
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcov -fprofile-arcs -ftest-coverage")

  # Link flags used for creating shared libraries
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lgcov -ftest-coverage")
endif()

# add custom directory to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules )

# store git hash
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# prepare configuration file
SET(PROGNAME "DEN2OBJ")
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "2")
SET(VERSION_MICRO "0")
message(STATUS "Writing configuration file in: ${CMAKE_CURRENT_SOURCE_DIR}/config.h")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/den2obj.pc.in ${CMAKE_BINARY_DIR}/den2obj.pc @ONLY)

# Enable release build
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# add OS specific
SET(BOOST_INCLUDEDIR "/usr/include")
SET(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")

find_package(OpenMP)
if (OPENMP_FOUND)
    option(HAS_OPENMP "OpenMP enabled" ON)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# set Boost
set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_STATIC_RUNTIME OFF)
set (BOOST_ALL_DYN_LINK OFF)

# Include libraries
find_package(PkgConfig REQUIRED)
find_package(Boost COMPONENTS regex iostreams filesystem REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(CPPUNIT REQUIRED) # for unit tests
pkg_check_modules(TCLAP tclap REQUIRED)
pkg_check_modules(EIGEN eigen3 REQUIRED)
pkg_check_modules(CRYPTO libcrypto REQUIRED) # for unit tests

# Set include folders
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${EIGEN_INCLUDE_DIRS}
                    ${TCLAP_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS}
                    ${CPPUNIT_INCLUDE_DIR})

# enable OpenVDB support
if(MOD_OPENVDB)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMOD_OPENVDB -I/opt/openvdb/include")
    SET(OPENVDB_LIBS "/opt/openvdb/lib/libopenvdb.a" "tbb" "Half" "blosc" "-lz")
    MESSAGE("Compiling with OpenVDB module")
endif()

# Add sources
file(GLOB SOURCES "*.cpp")
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/den2obj.cpp)
add_library(den2objsources STATIC ${SOURCES})
add_library(den2obj SHARED ${SOURCES})
add_executable(den2obj_exec ${CMAKE_CURRENT_SOURCE_DIR}/den2obj.cpp)
set_property(TARGET den2obj_exec PROPERTY OUTPUT_NAME den2obj)
file(GLOB DEN2OBJ_HEADERS "*.h")
set_target_properties(den2obj PROPERTIES PUBLIC_HEADER "${DEN2OBJ_HEADERS}")

# Set C++17
add_definitions(-std=c++17 -march=native)

# Link libraries
#SET(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=\$ORIGIN/lib")
target_link_libraries(den2obj_exec 
                      den2objsources 
                      ${OPENVDB_LIBS} 
                      ${Boost_LIBRARIES} 
                      ${LIBLZMA_LIBRARIES} 
                      ${ZLIB_LIBRARIES}
                      ${BZIP2_LIBRARIES})

if(NOT DISABLE_TEST)
    message("[USER] Performing unit tests")
    enable_testing ()
    add_subdirectory("test")
else()
    message("[USER] Testing routine disabled")
endif()

###
# Installing
##
include(GNUInstallDirs)

# install shared library and header files
install (TARGETS den2obj 
         DESTINATION bin
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/den2obj)

# install pkg-config files
install(FILES ${CMAKE_BINARY_DIR}/den2obj.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)