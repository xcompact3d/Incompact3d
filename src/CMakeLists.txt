include_directories(${CMAKE_SOURCE_DIR}/decomp2d
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
                    BC-ABL.f90
                    BC-Channel-flow.f90
                    BC-Cylinder.f90
                    BC-dbg-schemes.f90
                    BC-Jet.f90
                    BC-Lock-exchange.f90
                    BC-Mixing-layer.f90
                    BC-Periodic-hill.f90
		    BC-Sandbox.f90
                    BC-TBL.f90
                    BC-TGV.f90
                    BC-Uniform.f90
                    BC-User.f90
                    case.f90
                    constants.f90
                    derive.f90
                    dynstall.f90
                    dynstall_legacy.f90
                    filters.f90
                    forces.f90
                    genepsi3d.f90
                    ibm.f90
                    implicit.f90
                    les_models.f90
                    module_param.f90
                    navier.f90
                    parameters.f90
                    poisson.f90
                    probes.f90
                    schemes.f90
                    statistics.f90
                    time_integrators.f90
                    tools.f90
                    transeq.f90
                    turbine.f90
                    variables.f90
                    visu.f90
                    xcompact3d.f90)
#add_executable(xcompact3d ${files_xcompact})
target_include_directories(xcompact3d PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
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
