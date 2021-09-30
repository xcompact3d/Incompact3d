install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/TBL)
install(FILES input_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/TBL)
# If testing active add test for Cylinder case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TBL")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test.i3d DESTINATION ${case_dir})
endif()
