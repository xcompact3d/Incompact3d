macro(CreateMPITest test_dir
		 case
		 input_file
		 additional_files)
  message(STATUS "Add Verification Test (MPI run) ${case}")
  set(case_dir "${test_dir}/${case}")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY ${input_file} DESTINATION ${case_dir})
  set(local_list "")
  list(APPEND local_list ${additional_files})
  foreach(ff IN LISTS local_list)
    message(STATUS "${ff}")
    file(COPY ${ff} DESTINATION ${case_dir})
  endforeach()
  if(ADIOS2_FOUND)
    file(COPY adios2_config.xml DESTINATION ${case_dir})
  endif()
  add_test(NAME ${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:xcompact3d> ${input_file} WORKING_DIRECTORY ${case_dir})
endmacro()

macro(CreateMPIPythonTest 
		 test_dir
		 case
		 input_file
		 PyScript
		 ResuFile)
  # Test main param
  message(STATUS "Add Validation Test (MPI run + Python verification) ${case}")
  set(case_dir "${test_dir}/${case}")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY ${input_file} DESTINATION ${case_dir})
  if(ADIOS2_FOUND)
    file(COPY adios2_config.xml DESTINATION ${case_dir})
  endif()
  # Python
  file(COPY ${PyScript} DESTINATION ${case_dir})
  # Ref Solution
  set(RefFile "reference_${ResuFile}")
  file(COPY ${RefFile} DESTINATION ${case_dir})
  # set the MPI command
  set(cmd "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:xcompact3d> ${input_file}")
  add_test(NAME ${case}
           COMMAND ${CMAKE_COMMAND}
           -DCMD=${cmd}
	   -DPYSCRIPT=${PyScript}
	   -DRESUFILE=${ResuFile}
	   -DREFFILE=${RefFile}
	   -P ${CMAKE_SOURCE_DIR}/tests/cmake/X3D_runtests.cmake
           WORKING_DIRECTORY "${case_dir}")
endmacro()
