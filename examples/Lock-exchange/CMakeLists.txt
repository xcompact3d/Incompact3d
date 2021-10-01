install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Lock-exchange)
install(FILES input_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Lock-exchange)
# If testing active add test for Lock-exchange case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Lock-exchange")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test.i3d DESTINATION ${case_dir})
endif()
