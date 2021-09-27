# ABL Convective
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/ABL-Convective)
install(FILES input_convective_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/ABL-Convective)
# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/ABL-Convective")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_convective_test.i3d DESTINATION ${case_dir})
endif()
# ABL Neutral
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/ABL-Neutral)
install(FILES input_neutral_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/ABL-Neutral)
# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/ABL-Neutral")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_neutral_test.i3d DESTINATION ${case_dir})
endif()
# ABL Stable
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/ABL-Stable)
install(FILES input_stable_test.i3d DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/ABL-Stable)
# If testing active add test for TGV case
if (${BUILD_TESTING})
	set(case_dir "${test_dir}/ABL-Stable")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY input_stable_test.i3d DESTINATION ${case_dir})
endif()
