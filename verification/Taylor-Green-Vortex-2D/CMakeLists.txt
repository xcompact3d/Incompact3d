include_directories(${PROJECT_BINARY_DIR}/decomp2d
                    ${PROJECT_BINARY_DIR}/src
                    ${PROJECT_BINARY_DIR}/verification/lib)

add_pfunit_ctest (tgv2D_00
  TEST_SOURCES init00.pf test00.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (tgv2D_00_imp1
  TEST_SOURCES init00_imp1.pf test00_imp1.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (tgv2D_11
  TEST_SOURCES init11.pf test11.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (tgv2D_11_imp2
  TEST_SOURCES init11_imp2.pf test11_imp2.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)
