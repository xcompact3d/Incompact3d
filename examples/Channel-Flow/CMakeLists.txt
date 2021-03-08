install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel-Flow)
install(FILES input.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel-Flow)
# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Channel-Flow")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input.i3d DESTINATION ${case_dir})
endif()
