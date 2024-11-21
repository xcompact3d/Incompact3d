set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(X3D_CREATE_TEST)

set(test_dir "${PROJECT_BINARY_DIR}/RunTests")

# Python 
find_package(Python3 REQUIRED)
if(NOT Python3_FOUND)
  message(FATAL_ERROR "Python3 support NOT FOUND! Python3 is needed to compare the solutions for the TGV case.")
endif()

# We need to find a way to check if 2DECOMP&FFT is build with ADIOS2
option(USE_ADIOS2 "Build XCompact with ADIOS2 library" OFF)

# TGV is the only test active by default
add_subdirectory(TGV-Taylor-Green-vortex)
if (${BUILD_TESTING_FULL})
  add_subdirectory(ABL-Atmospheric-Boundary-Layer)
  add_subdirectory(Cavity)
  add_subdirectory(Channel)
  if(NOT ${USE_ADIOS2})
    add_subdirectory(Cylinder-wake)
  endif()
  add_subdirectory(Gravity-current)
  add_subdirectory(MHD)
  add_subdirectory(Mixing-layer)
  add_subdirectory(Particle-Tracking)
  add_subdirectory(Periodic-hill)
  add_subdirectory(TBL-Turbulent-Boundary-Layer)
  add_subdirectory(Wind-Turbine)
endif (${BUILD_TESTING_FULL})


