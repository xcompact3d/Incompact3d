install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
# if testing active create a working dir for testing
if (${BUILD_TESTING})
  file(MAKE_DIRECTORY ${test_dir})
endif()
add_subdirectory(Taylor-Green-Vortex)
add_subdirectory(Channel-Flow)
