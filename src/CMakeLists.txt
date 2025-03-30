include_directories(${2decomp_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(xcompact3d 
                    acl_controller.f90 
                    acl_elem.f90 
                    acl_farm_controller.f90 
                    acl_model.f90
                    acl_out.f90  
                    acl_source.f90
                    acl_turb.f90
                    acl_utils.f90 
                    adm.f90
                    airfoils.f90
                    BC-Ellipsoid.f90
                    Case-ABL.f90
                    Case-Cavity.f90
                    Case-Channel.f90
                    Case-Cylinder-wake.f90
                    Case-Generic.f90
                    Case-Gravity-current.f90
                    Case-Mixing-layer.f90
                    Case-Periodic-hill.f90
                    Case-Pipe-flow.f90
                    Case-Sandbox.f90
                    Case-TBL.f90
                    Case-TGV.f90
                    Case-Uniform.f90
                    Case-PTBL.f90
                    case.f90
                    derive.f90
                    dynstall.f90
                    dynstall_legacy.f90
                    ellip_utils.f90
                    filters.f90
                    forces.f90
                    genepsi3d.f90
                    ibm.f90
                    implicit.f90
                    les_models.f90
                    mhd.f90
                    module_param.f90
                    mptool.f90
                    navier.f90
                    parameters.f90
                    particle.f90
                    poisson.f90
                    probes.f90
                    schemes.f90
                    statistics.f90
		    stretching.f90
                    time_integrators.f90
                    tools.f90
                    transeq.f90
                    turbine.f90
                    utilities.f90
                    variables.f90
                    visu.f90
                    xcompact3d.f90)

target_link_libraries(xcompact3d PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(ADIOS2_FOUND)
  #target_link_libraries(xcompact3d ${ADIOS_LIBRARIES})
  target_link_libraries(xcompact3d PRIVATE adios2::fortran_mpi)
endif()

install(TARGETS xcompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
