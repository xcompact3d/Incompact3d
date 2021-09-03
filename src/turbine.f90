module turbine

contains

  subroutine init_turbines(ux1,uy1,uz1)

    USE decomp_2d      
    USE actuator_line_model
    USE actuator_line_source
    USE actuator_disc_model
    USE param

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    ! Initialise the turbine models
    if (iturbine == 1) then
      call actuator_line_model_init(Nturbines,Nactuatorlines,TurbinesPath,ActuatorlinesPath,dt)
      call initialize_actuator_source
      call Compute_Momentum_Source_Term_pointwise
      if (nrank==0.and.mod(itime,iturboutput)==0) then
        call actuator_line_model_write_output(itime/iturboutput)
      endif
    else if (iturbine == 2) then
      call actuator_disc_model_init(Ndiscs,admCoords,C_T,aind)
      call actuator_disc_model_compute_source(ux1,uy1,uz1)
    endif

  end subroutine init_turbines

  subroutine compute_turbines(ux1,uy1,uz1)

    USE decomp_2d      
    USE actuator_line_model
    USE actuator_line_source
    USE actuator_disc_model
    USE param
   
    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
     
    if (iturbine == 1) then
      if (nrank==0) then
        write(6,*) ''
        write(6,*) 'Unsteady Actuator Line Model INFO:'
      endif
      call Compute_Momentum_Source_Term_pointwise
      call actuator_line_model_update(t,dt)
      if (nrank==0) then
        write(6,*) ''
      endif
    else if (iturbine == 2) then
      call actuator_disc_model_compute_source(ux1,uy1,uz1)
    endif

  end subroutine compute_turbines

  subroutine turbine_output()

    USE param
    USE actuator_line_model
    USE actuator_disc_model
    USE decomp_2d_io

    implicit none   

    if (nrank==0.and.mod(itime,iturboutput)==0) then
      if (iturbine == 1) then
        call actuator_line_model_write_output(itime/iturboutput)
      else if (iturbine == 2) then
        call actuator_disc_model_write_output(itime/iturboutput)
      endif
    endif

  end subroutine turbine_output

end module turbine
