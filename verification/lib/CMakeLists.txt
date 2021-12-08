include_directories(${PROJECT_BINARY_DIR}/decomp2d
                    ${PROJECT_BINARY_DIR}/src)

set(SRCFILES verif_TGV2D.f90 x3d_unit_testing_tools.f90)

add_library(x3dverif STATIC ${SRCFILES})

add_dependencies(x3dverif decomp2d x3d)
target_link_libraries(x3dverif PRIVATE decomp2d x3d)
if (MPI_FOUND)
	target_link_libraries(x3dverif PRIVATE MPI::MPI_Fortran)
endif(MPI_FOUND)
if (ADIOS2_FOUND)
	target_link_libraries(x3dverif PRIVATE adios2::fortran_mpi)
endif(ADIOS2_FOUND)
