cmake_minimum_required(VERSION 3.20)                        
cmake_policy(SET CMP0074 NEW)

project(Xcompact3d LANGUAGES Fortran)
set(version 5.0)
if (IO_BACKEND MATCHES "adios2")
  # Can be useful to also activate CXX, sometimes is needed by packages
  enable_language(C CXX)
endif (IO_BACKEND MATCHES "adios2")

set(DESCRIPTION "Building Xcompact3d using cmake")

message(STATUS "building ${PROJECT_NAME}")

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/opt" CACHE PATH "..." FORCE)
endif()

# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
       "Choose the type of build, options are: NONE DEV DEBUG RELEASE."
       FORCE)
endif (NOT CMAKE_BUILD_TYPE)
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
if (CMAKE_BUILD_TYPE MATCHES "NONE")
  message (STATUS "Selected build type : None")
elseif (CMAKE_BUILD_TYPE MATCHES "RELEASE")
  message (STATUS "Selected build type : Release")
elseif (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  message (STATUS "Selected build type : Debug")
elseif (CMAKE_BUILD_TYPE MATCHES "DEV")
  message (FATAL_ERROR "DEV build is a work in progress and currently not supported by Xcompact3d")
else (CMAKE_BUILD_TYPE MATCHES "NONE")
  message (FATAL_ERROR "Invalid build type, options are: NONE DEV DEBUG RELEASE.")
endif (CMAKE_BUILD_TYPE MATCHES "NONE")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/compilers" "${CMAKE_SOURCE_DIR}/cmake/fft")
find_package(DECOMP2D REQUIRED)

# Despite IO is handled by 2DECOMP&FFT we need to set it
set (IO_BACKEND "mpi" CACHE STRING "Default IO backend (mpi (default) or adios2)")
set_property(CACHE IO_BACKEND PROPERTY STRINGS mpi adios2)

include(X3D_MPI)

include(X3D_Compilers)

# Add the sources
add_subdirectory(src)

# Add tests
option(BUILD_TESTING "Build with test and use only TGV case" ON)
option(BUILD_TESTING_PARTICLE "Build with test and use only two cases with particle tracking" ON)
option(BUILD_TESTING_FULL "Build all tests" OFF)
if (${BUILD_TESTING})
  include(CTest)
  add_subdirectory(examples)
  add_subdirectory(tests) 
endif (${BUILD_TESTING})

# Add a prettify target
#add_custom_target(format sh ${CMAKE_SOURCE_DIR}/scripts/format.sh
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})