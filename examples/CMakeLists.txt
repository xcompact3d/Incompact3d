install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

# if testing active create a working dir for testing
set(test_dir "${PROJECT_BINARY_DIR}/RunTests")
file(MAKE_DIRECTORY ${test_dir})

# We need to find a way to check if 2DECOMP&FFT is build with ADIOS2
option(USE_ADIOS2 "Build XCompact with ADIOS2 library" OFF)

option(BUILD_TESTING_FULL "Build with tests" OFF)
add_subdirectory(TGV)

if (${BUILD_TESTING_FULL})
  add_subdirectory(ABL)
  add_subdirectory(ABL-old)
  add_subdirectory(Cavity)
  add_subdirectory(Channel-Flow)
  if(NOT ${USE_ADIOS2})
    add_subdirectory(Cylinder)
  endif()
  add_subdirectory(Lock-exchange)
  add_subdirectory(Mixing-layer)
  add_subdirectory(Turbulent-Boundary-Layer)
  add_subdirectory(Wind-Turbine/NREL-5MW_ALM)
endif (${BUILD_TESTING_FULL})
