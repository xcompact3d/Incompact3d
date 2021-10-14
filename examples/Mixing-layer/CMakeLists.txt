install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Mixing-layer)
install(FILES input_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Mixing-layer)
# If testing active add test for Mixing-layer case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Mixing-layer")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test.i3d DESTINATION ${case_dir})
endif()
