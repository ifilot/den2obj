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
cmake_minimum_required(VERSION 2.8.12)
project (den2obj)

# add custom directory to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules )

# prepare configuration file
SET(PROGNAME "DEN2OBJ")
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "6")
SET(VERSION_MICRO "0")
configure_file(config.h.in ../src/config.h @ONLY)

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
find_package(CPPUNIT REQUIRED) # for unit tests
pkg_check_modules(TCLAP tclap REQUIRED)
pkg_check_modules(EIGEN eigen3 REQUIRED)

# Set include folders
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${EIGEN_INCLUDE_DIRS}
                    ${TCLAP_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIRS}
                    ${CPPUNIT_INCLUDE_DIR})

# Add sources
file(GLOB SOURCES "*.cpp")
add_executable(den2obj ${SOURCES})
add_library(den2objsources STATIC ${SOURCES})

# Set C++17
add_definitions(-std=c++17 -march=native)

# Link libraries
#SET(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=\$ORIGIN/lib")
target_link_libraries(den2obj ${Boost_LIBRARIES} openvdb tbb Half)

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
install (TARGETS den2obj DESTINATION bin)
