include_directories(${PROJECT_BINARY_DIR}/decomp2d
                    ${PROJECT_BINARY_DIR}/src
                    ${PROJECT_BINARY_DIR}/verification/lib)

add_pfunit_ctest (verif_x_00
  TEST_SOURCES initx00.pf testx00.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_00
  TEST_SOURCES inity00.pf testy00.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_00_str
  TEST_SOURCES inity00_str.pf testy00_str.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_z_00
  TEST_SOURCES initz00.pf testz00.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_x_11
  TEST_SOURCES initx11.pf testx11.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_11
  TEST_SOURCES inity11.pf testy11.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_11_str
  TEST_SOURCES inity11_str.pf testy11_str.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_z_11
  TEST_SOURCES initz11.pf testz11.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_x_22
  TEST_SOURCES initx22.pf testx22.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_22
  TEST_SOURCES inity22.pf testy22.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_22_str
  TEST_SOURCES inity22_str.pf testy22_str.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_z_22
  TEST_SOURCES initz22.pf testz22.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_x_12
  TEST_SOURCES initx12.pf testx12.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_12
  TEST_SOURCES inity12.pf testy12.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_12_str
  TEST_SOURCES inity12_str.pf testy12_str.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_z_12
  TEST_SOURCES initz12.pf testz12.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_x_21
  TEST_SOURCES initx21.pf testx21.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_21
  TEST_SOURCES inity21.pf testy21.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_y_21_str
  TEST_SOURCES inity21_str.pf testy21_str.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)

add_pfunit_ctest (verif_z_21
  TEST_SOURCES initz21.pf testz21.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)
