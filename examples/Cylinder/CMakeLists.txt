install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Cylinder)
install(FILES input_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Cylinder)
# If testing active add test for Cylinder case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Cylinder")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test.i3d DESTINATION ${case_dir})
endif()
