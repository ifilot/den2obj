# Examples

This folder contains a number of scripts how to use Den2Obj as a (shared)
library in your own programs and applications.

## Requirements

To use Den2Obj as a shared library in your applications, you need to dynamically
link against `libden2obj.so` as well as against a number of required libraries:

* Boost (regex, iostreams and filesystem)
* GZIP
* LZMA
* BZ2

Besides these libraries, there is also a header-only dependency on the Eigen3
library.

A convenient way to ensure this in your application is by making use of
`CMake` and using the following settings in your `CmakeLists.txt file`

```cmake
# use Boost
SET(BOOST_INCLUDEDIR "/usr/include")
SET(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")

# set Boost
set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
set (BOOST_ALL_DYN_LINK OFF)

# use OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    option(HAS_OPENMP "OpenMP enabled" ON)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# include libraries
find_package(PkgConfig REQUIRED)
find_package(Boost COMPONENTS regex iostreams filesystem REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
pkg_check_modules(DEN2OBJ den2obj REQUIRED)
pkg_check_modules(EIGEN eigen3 REQUIRED)
```

and finally add the required libraries to your executable

```cmake
target_link_libraries(den2obj-shared-example 
                      ${DEN2OBJ_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${LIBLZMA_LIBRARIES} 
                      ${ZLIB_LIBRARIES}
                      ${BZIP2_LIBRARIES})
```

An example of this is provided in [shared/CMakeLists.txt](shared/CMakeLists.txt).

## Compilation instructions

```bash
mkdir build && cd build
cmake ../shared
make -j
```