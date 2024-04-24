install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

# if testing active create a working dir for testing
set(test_dir "${PROJECT_BINARY_DIR}/RunTests")
file(MAKE_DIRECTORY ${test_dir})


if (${BUILD_TESTING_FULL})
  add_subdirectory(ABL-Atmospheric-Boundary-Layer)
  add_subdirectory(Cavity)
  add_subdirectory(Channel)
  if(NOT ${USE_ADIOS2})
    add_subdirectory(Cylinder-wake)
  endif()
  add_subdirectory(Gravity-current)
  add_subdirectory(Mixing-layer)
  add_subdirectory(OTV-Ozage-Tang-vortex-2D)
  add_subdirectory(TBL-Turbulent-Boundary-Layer)
  add_subdirectory(Wind-Turbine/NREL-5MW_ALM)
endif (${BUILD_TESTING_FULL})
