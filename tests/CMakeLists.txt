set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(X3D_CREATE_TEST)

set(test_dir "${PROJECT_BINARY_DIR}/RunTests")

# TGV is from here
# TGV is the only test active by default
add_subdirectory(TGV-Taylor-Green-vortex)

