# As part of installation we'll copy the examples
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(X3D_ADD_EXAMPLE)
# Tests are in the folder tests
# This is only to have in the installation director the examples
set(examples_dir "${CMAKE_INSTALL_PREFIX}/examples")
install(DIRECTORY DESTINATION ${example_dir})

add_subdirectory(ABL-Atmospheric-Boundary-Layer)
add_subdirectory(Cavity)
add_subdirectory(Channel)
add_subdirectory(Cylinder-wake)
add_subdirectory(Gravity-current)
add_subdirectory(MHD)
add_subdirectory(Mixing-layer)
add_subdirectory(Periodic-hill)
add_subdirectory(Pipe-Flow)
add_subdirectory(TBL-Turbulent-Boundary-Layer)
add_subdirectory(TGV-Taylor-Green-vortex)
add_subdirectory(Wind-Turbine)
