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

cmake_minimum_required(VERSION 3.16)
project(den2obj VERSION 1.3.0 LANGUAGES CXX)

option(USE_GCOV "Build with gcov coverage instrumentation" OFF)
option(MOD_OPENVDB "Enable OpenVDB support" OFF)
option(HAS_OPENMP "Enable OpenMP support when available" ON)
option(DISABLE_TEST "Disable building tests" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(PROGNAME "DEN2OBJ")
set(VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(VERSION_MICRO "${PROJECT_VERSION_PATCH}")

execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

message(STATUS "Writing configuration file in: ${CMAKE_CURRENT_BINARY_DIR}/config.h")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/den2obj.pc.in ${CMAKE_CURRENT_BINARY_DIR}/den2obj.pc @ONLY)

find_package(PkgConfig REQUIRED)
find_package(Boost REQUIRED COMPONENTS regex iostreams filesystem)
find_package(LibLZMA REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
pkg_check_modules(TCLAP REQUIRED IMPORTED_TARGET tclap)
pkg_check_modules(EIGEN REQUIRED IMPORTED_TARGET eigen3)

if(HAS_OPENMP)
    find_package(OpenMP)
endif()

set(OPENVDB_LIBS "")
if(MOD_OPENVDB)
    set(OPENVDB_LIBS openvdb tbb blosc ZLIB::ZLIB)
    message(STATUS "Compiling with OpenVDB module")
endif()

file(GLOB SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/den2obj.cpp")
file(GLOB DEN2OBJ_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(den2objsources STATIC ${SOURCES})
add_library(den2obj SHARED ${SOURCES})
add_executable(den2obj_exec ${CMAKE_CURRENT_SOURCE_DIR}/den2obj.cpp)
set_property(TARGET den2obj_exec PROPERTY OUTPUT_NAME den2obj)
set_target_properties(den2obj PROPERTIES PUBLIC_HEADER "${DEN2OBJ_HEADERS}")

foreach(target den2objsources den2obj den2obj_exec)
    target_compile_features(${target} PUBLIC cxx_std_17)
    target_compile_options(${target} PRIVATE -march=native)
    target_include_directories(${target}
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}
            ${Boost_INCLUDE_DIRS}
            ${TCLAP_INCLUDE_DIRS}
            ${EIGEN_INCLUDE_DIRS}
    )

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
    endif()

    if(MOD_OPENVDB)
        target_compile_definitions(${target} PUBLIC MOD_OPENVDB)
    endif()
endforeach()

if(USE_GCOV)
    foreach(target den2objsources den2obj den2obj_exec)
        target_compile_options(${target} PRIVATE --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endforeach()
endif()

target_link_libraries(den2objsources
    PUBLIC
        Boost::regex
        Boost::iostreams
        Boost::filesystem
        LibLZMA::LibLZMA
        ZLIB::ZLIB
        BZip2::BZip2
        ${OPENVDB_LIBS}
)

target_link_libraries(den2obj PUBLIC den2objsources)
target_link_libraries(den2obj_exec PRIVATE den2objsources)

include(CTest)
if(DISABLE_TEST)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
endif()

if(BUILD_TESTING)
    message(STATUS "[USER] Performing unit tests")
    add_subdirectory(test)
else()
    message(STATUS "[USER] Testing routine disabled")
endif()

include(GNUInstallDirs)

install(TARGETS den2obj den2obj_exec
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/den2obj)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/den2obj.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
