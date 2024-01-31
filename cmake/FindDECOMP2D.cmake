# - Find the 2decomp-fft library
find_package(decomp2d
             PATHS ${CMAKE_SOURCE_DIR}/decomp2d/build)
if (decomp2d_FOUND)
  message(STATUS "2decomp-fft FOUND")
else(decomp2d_FOUND)
  message(STATUS "2decomp-fft PATH not available we'll try to download and install")
  configure_file(${CMAKE_SOURCE_DIR}/cmake/decomp2d/downloadBuild2decomp.cmake.in decomp2d-build/CMakeLists.txt)
  #message("Second CMAKE_GENERATOR ${CMAKE_GENERATOR}") 
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build )
  if(result)
      message(FATAL_ERROR "CMake step for 2decomp-fft failed: ${result}")
  else()
      message("CMake step for 2decomp-fft completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
         RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build )
  if(result)
      message(FATAL_ERROR "Build step for 2decomp-fft failed: ${result}")
  endif()
  set(D2D_ROOT ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build/downloadBuild2decomp-prefix/src/downloadBuild2decomp-build)
  find_package(decomp2d REQUIRED
	  PATHS ${D2D_ROOT})
endif(decomp2d_FOUND)



