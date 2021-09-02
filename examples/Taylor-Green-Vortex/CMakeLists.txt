install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor-Green-Vortex)
install(FILES input_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor-Green-Vortex)
# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TGV")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test.i3d DESTINATION ${case_dir})
endif()
