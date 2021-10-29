subroutine GebraadController(Turbine,Ntur,time)
    
    use actuator_line_model_utils
    use actuator_line_turbine
    use actuator_line_controller
    use param, only: zero, one
    use constants

    implicit none
    real(mytype) :: time
    integer :: Ntur
    type(TurbineType), dimension(Ntur), intent(inout) :: turbine
    integer :: i,j,k
    real(mytype) :: TotalPower, WSRotorAve

    TotalPower=zero
    do i=1,Ntur
       TotalPower=TotalPower+Turbine(i)%Power
       call compute_rotor_upstream_velocity(Turbine(i))
    enddo

    do i=1,Ntur
       Turbine(i)%yaw_angle=zero
       Turbine(i)%shaft_tilt_angle=zero

       ! APPLY FARM LEVEL CONTROL

       ! Rotate the turbine according to the tilt and yaw angle
       ! Yaw
       call rotate_turbine(turbine(i),(/zero,one,zero/),turbine(i)%yaw_angle*conrad)
       ! Tilt
       call rotate_turbine(turbine(i),(/zero,zero,one/),-turbine(i)%shaft_tilt_angle*conrad)
   
       ! Set the rotational axis
       call QuatRot(turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3),turbine(i)%yaw_angle*conrad,zero,one,zero,zero,zero,zero,&
                    turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3))
       call QuatRot(turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3),-turbine(i)%shaft_tilt_angle*conrad,zero,zero,one,zero,zero,zero,&
                    turbine(i)%RotN(1),turbine(i)%RotN(2),turbine(i)%RotN(3))
    enddo

end subroutine GebraadController        
