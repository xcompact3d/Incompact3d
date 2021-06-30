install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel-Flow-X)
install(FILES input_test_x.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel-Flow-X)
# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Channel-Flow-X")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test_x.i3d DESTINATION ${case_dir})
endif()
# Z dir test
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel-Flow-Z)
install(FILES input_test_z.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel-Flow-Z)
# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Channel-Flow-Z")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_test_z.i3d DESTINATION ${case_dir})
endif()
