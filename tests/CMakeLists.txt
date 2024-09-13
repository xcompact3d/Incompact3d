set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(X3D_CREATE_TEST)

set(test_dir "${PROJECT_BINARY_DIR}/RunTests")

# Python 
find_package(Python3 REQUIRED)
if(NOT Python3_FOUND)
  message(FATAL_ERROR "Python3 support NOT FOUND! Python3 is needed to compare the solutions for the TGV case.")
endif()

# TGV is from here
# TGV is the only test active by default
add_subdirectory(TGV-Taylor-Green-vortex)

