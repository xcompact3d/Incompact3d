# TGV
set(case "TGV-Taylor-Green-vortex")
set(case_dir "${test_dir}/${case}")
file(MAKE_DIRECTORY ${case_dir})
set(input_file "reference_input.i3d")
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/${case})
install(FILES ${input_file} DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/${case})
install(FILES reference_time_evol.dat DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/${case})
install(FILES ../../compare_TGV_time_evolution.py DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/${case})
# If testing active add test for TGV case
file(COPY ${input_file} DESTINATION ${case_dir})
file(COPY reference_time_evol.dat DESTINATION ${case_dir})
file(COPY ../../compare_TGV_time_evolution.py DESTINATION ${case_dir})
add_test(NAME ${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:xcompact3d> ${input_file} WORKING_DIRECTORY ${case_dir})
