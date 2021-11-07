cmake_minimum_required(VERSION 3.0.2)                        
cmake_policy(SET CMP0048 NEW)

project(xcompact3d LANGUAGES Fortran)

set(AUTHOR "Stefano Rolfo;Charles Moulinec")
set(AUTHOR_DETAILS "stefano.rolfo@stfc.ac.uk;charles.moulinec@stfc.ac.uk")
set(DESCRIPTION "Building xcompact3d using cmake")

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

# Find the modules included with Xcompact
#SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} )
message(STATUS "COMP ID ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler name ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler version ${CMAKE_Fortran_COMPILER_VERSION}")
if (Fortran_COMPILER_NAME MATCHES "GNU")
  # gfortran
  message(STATUS "Setting gfortran flags")
  set(CMAKE_Fortran_FLAGS "-cpp -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none")
  if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
    message(STATUS "Set New Fortran basic flags")
    set(CMAKE_Fortran_FLAGS "-cpp -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none -fallow-argument-mismatch")
  endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-DDEBG -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "Intel")
  message(STATUS "Setting ifort flags")
  set(CMAKE_Fortran_FLAGS "-fpp -xHost -heaparrays -safe-cray-ptr -g -traceback")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -DDEBG")
  #set(CMAKE_Fortran_FLAGS "-cpp xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr")
elseif (Fortran_COMPILER_NAME MATCHES "NAG")
  message(STATUS "Setting nagfor flags")
  set(CMAKE_Fortran_FLAGS "-fpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "Cray")
  message(STATUS "Setting cray fortran flags")
  set(CMAKE_Fortran_FLAGS "-eF -g -N 1023")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "PGI")
  message(STATUS "Setting PGI fortran flags")
  set(CMAKE_Fortran_FLAGS "-cpp -acc -Mfree -Kieee -Minfo=accel -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -DDEBG")
elseif (Fortran_COMPILER_NAME MATCHES "Fujitsu")
  message(STATUS "Setting Fujitsu fortran flags")
  set (CMAKE_Fortran_FLAGS "-Cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
else (Fortran_COMPILER_NAME MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "GNU")

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions("-DDEBG")
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")


find_package(MPI REQUIRED)
# Stop if there is no MPI_Fortran_Compiler
if (MPI_Fortran_COMPILER)
    message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
else (MPI_Fortran_COMPILER)
    message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)
# Warning if Include are not found => can be fixed with more recent cmake version
if (MPI_FOUND)
    message(STATUS "MPI FOUND: ${MPI_FOUND}")
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
else (MPI_FOUND)
    message(STATUS "NO MPI include have been found. The executable won't be targeted with MPI include")
    message(STATUS "Code will compile but performaces can be compromised")
    message(STATUS "Using a CMake vers > 3.10 should solve the problem")
    message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (MPI_FOUND)

execute_process(
  COMMAND git describe --tag --long --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DVERSION=\"${GIT_VERSION}\"")
option(DOUBLE_PRECISION "Build Xcompact with double precision" ON)
if (DOUBLE_PRECISION)
  add_definitions("-DDOUBLE_PREC")
endif()

option(SINGLE_PRECISION_OUTPUT "Build XCompact with output in single precision" OFF)
if (SINGLE_PRECISION_OUTPUT)
  add_definitions("-DSAVE_SINGLE")
endif()

# FFT options
set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw3 mkl)

# Look for ADIOS2 library
option(USE_ADIOS2 "Build XCompact with ADIOS2 library" OFF)
if(USE_ADIOS2)
  #find_package(ADIOS2 COMPONENTS fortran)
  find_package(ADIOS2)
endif()
if(ADIOS2_FOUND)
  message(STATUS "ADIOS INCLUDE ${ADIOS2_INCLUDE_DIRS}")
  include_directories(${ADIOS2_INCLUDE_DIRS})
  add_definitions("-DADIOS2")
else(ADIOS_FOUND)
  message(STATUS "ADIOS2 not found")
endif()

# Create a static library for the fft
add_subdirectory(decomp2d)

# Create the Xcompact3d executable
add_subdirectory(src)

# Create the Xcompact3d executable
add_subdirectory(post_vtk)

# Create an example dir with all input.i3d example files
option(BUILD_TESTING "Build with tests" ON)
set(test_dir "${PROJECT_BINARY_DIR}/Test")
add_subdirectory(examples)
message(STATUS "Before test main ${test_dir}")
if (${BUILD_TESTING})
  include(CTest)
  message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
  message(STATUS "MPI EXEC: ${MPIEXEC_EXECUTABLE}")
  add_test(NAME TGvortex COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test.i3d WORKING_DIRECTORY ${test_dir}/TGV)
  add_test(NAME Channel-Flow-X COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test_x.i3d WORKING_DIRECTORY ${test_dir}/Channel-Flow-X)
  add_test(NAME Channel-Flow-Z COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test_z.i3d WORKING_DIRECTORY ${test_dir}/Channel-Flow-Z)
  if(NOT ${USE_ADIOS2})
    add_test(NAME Cylinder COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test.i3d WORKING_DIRECTORY ${test_dir}/Cylinder)
    add_test(NAME CylinderMoving COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test_moving.i3d WORKING_DIRECTORY ${test_dir}/CylinderMoving)
  endif()
  add_test(NAME Lock-exchange COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test.i3d WORKING_DIRECTORY ${test_dir}/Lock-exchange)
  add_test(NAME Mixing-layer COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test.i3d WORKING_DIRECTORY ${test_dir}/Mixing-layer)
  add_test(NAME ABL-Convective COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_convective_test.i3d WORKING_DIRECTORY ${test_dir}/ABL-Convective)
  add_test(NAME ABL-Neutral COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_neutral_test.i3d WORKING_DIRECTORY ${test_dir}/ABL-Neutral)
  add_test(NAME ABL-Stable COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_stable_test.i3d WORKING_DIRECTORY ${test_dir}/ABL-Stable)
  add_test(NAME TBL COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test.i3d WORKING_DIRECTORY ${test_dir}/TBL)
  add_test(NAME Wind-Turbine COMMAND ${MPIEXEC_EXECUTABLE} -n ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_INSTALL_PREFIX}/bin/xcompact3d input_test.i3d WORKING_DIRECTORY ${test_dir}/Wind-Turbine)
endif()

