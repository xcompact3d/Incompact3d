include_directories(${CMAKE_SOURCE_DIR}/post_vtk)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(xcompact3d_paraview_vtk
	       paraview_incompact3d_vtk.f90)

install(TARGETS xcompact3d_paraview_vtk
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
