install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)

# if testing active create a working dir for testing
set(test_dir "${PROJECT_BINARY_DIR}/RunTests")
file(MAKE_DIRECTORY ${test_dir})
# TGV is from here
add_subdirectory(data/Taylor-Green-Vortex)

