subroutine GebraadController(Turbine,Ntur,time)
    
    use actuator_line_model_utils
    use actuator_line_turbine
    use actuator_line_controller

    implicit none
    real(mytype) :: time
    integer :: Ntur
    type(TurbineType), dimension(Ntur), intent(inout) :: turbine
    integer :: i,j,k
    real(mytype) :: TotalPower, WSRotorAve

    TotalPower=0.0_mytype
    do i=1,Ntur
    TotalPower=TotalPower+Turbine(i)%Power
    call compute_rotor_upstream_velocity(Turbine(i))
    enddo

    do i=1,Ntur
    Turbine(i)%yaw_angle=0.0_mytype
    Turbine(i)%shaft_tilt_angle=0.0_mytype

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! APPLY FARM LEVEL CONTROL
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ! Rotate the turbine according to the tilt and yaw angle
    ! Yaw
    call rotate_turbine(turbine(i),(/0.0d0,1.0d0,0.0d0/),turbine(i)%yaw_angle*pi/180.0d0)
    ! Tilt
    call rotate_turbine(turbine(i),(/0.0d0,0.0d0,1.0d0/),-turbine(i)%shaft_tilt_angle*pi/180.0d0)
   
    ! Set the rotational axis
    call QuatRot(turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3),turbine(i)%yaw_angle*pi/180.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.d0,&
            turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3))
    call QuatRot(turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3),-turbine(i)%shaft_tilt_angle*pi/180.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.d0,&
            turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3))
    enddo

end subroutine GebraadController        
