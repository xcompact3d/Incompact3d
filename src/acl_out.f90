module actuator_line_write_output
    
    use actuator_line_turbine
    use actuator_line_element
    use actuator_line_controller

contains
 
    !*******************************************************************************
    !
    subroutine actuator_line_element_write_output(act_line,dir)
    !
    !*******************************************************************************

      implicit none
      type(ActuatorLineType), intent(in) :: act_line
      character(len=100), intent(in) :: dir
      character(len=22) :: Format
      integer :: ielem

      open(2017,File=trim(dir)//'_'//trim(act_line%name)//'.load')
      write(2017,*) 'ielem,X,Y,Z,rdist/R,AOA,adot,RE,ur,CL,CD,CM25,Cn,Ct,Fn,Ft,F1'
      write(2017,*) '[-], [m],[m],[m],[-],[deg],[deg/s],[-],[m/s],[-],[-],[-],[-],[-],[N],[N],[-]'
      Format="(I5,A,19(E14.7,A))"
      do ielem=1,act_line%NElem
         write(2017,Format) ielem,',',act_line%PEx(ielem),',',act_line%PEy(ielem),',',act_line%PEz(ielem),',',act_line%ERdist(ielem)/act_line%L,',',act_line%EAOA(ielem)*180/pi,',',act_line%EAOAdot(ielem)*pi/180,',',act_line%ERE(ielem),',',act_line%EUr(ielem),',',act_line%ECL(ielem),',',act_line%ECD(ielem),',',act_line%ECM(ielem),',',act_line%ECN(ielem),',',act_line%ECT(ielem),',',act_line%EFn(ielem),',',act_line%EFt(ielem),',',act_line%EEndeffects_factor(ielem)
      enddo
      close(2017)
 
    end subroutine actuator_line_element_write_output
    
    !*******************************************************************************
    !
    subroutine dynamic_stall_write_output(act_line,dir)
    !
    !*******************************************************************************

      implicit none
      type(ActuatorLineType), intent(in) :: act_line
      character(len=100), intent(in) :: dir
      character(len=22) :: Format
      integer :: ielem

      open(2018,File=trim(dir)//'_'//trim(act_line%name)//'.dynstall')
      write(2018,*) 'ielem,rdist/R,pitch,AOA,f'
      Format="(I5,A,15(E14.7,A))"
      do ielem=1,act_line%NElem
         write(2018,Format)ielem,',',act_line%ERdist(ielem)/act_line%L,',',act_line%Epitch(ielem)*180/pi,',',act_line%EAOA(ielem)*180/pi,',',act_line%EDynstall(ielem)%fprime
      enddo
      close(2018)
 
    end subroutine dynamic_stall_write_output
    
    !*******************************************************************************
    !
    subroutine actuator_line_turbine_write_output(turbine,dir)
    !
    !*******************************************************************************

      implicit none
      type(TurbineType), intent(in) :: turbine
      character(len=100), intent(in) :: dir
      character(len=22) :: Format
        
      open(2016,File=trim(dir)//'_'//trim(turbine%name)//'.perf')
      write(2016,*) 'Number of Revs, GeneratorSpeed, GeneratorTorque, BladePitch1, BladePitch2, BladePitch3, Omega, DOmega, Ux, Uy, Uz, Thrust, Torque, Power'
      write(2016,*) '[-], [rad/s], [N m], [deg], [deg], [deg], [rad/s], [rad/s], [m/s], [m/s], [m/s], [N], [N m], [W]'
      Format="(14(E14.7,A))" 
      write(2016,Format) turbine%AzimAngle/(2*pi),',',turbine%controller%GenSpeed,',',turbine%controller%GenTrq,',',turbine%controller%PitCom(1)*180.0_mytype/pi,',',turbine%controller%PitCom(2)*180.0_mytype/pi,',',turbine%controller%PitCom(3)*180.0_mytype/pi,',',turbine%angularVel,',',turbine%deltaOmega,',',turbine%Ux_upstream,',',turbine%Uy_upstream,',',turbine%Uz_upstream,',',turbine%Thrust,',',turbine%Torque,',',turbine%Power
      close(2016)

    end subroutine actuator_line_turbine_write_output   
    
    !*******************************************************************************
    !
    subroutine actuator_line_turbine_write_statistics(turbine,dir)
    !
    !*******************************************************************************

      implicit none
      type(TurbineType), intent(in) :: turbine
      character(len=100), intent(in) :: dir
        
      open(2019,File=trim(dir)//'_'//trim(turbine%name)//'.stat')
      write(2019,*) 'T_ave, P_ave, Torque_ave'
      write(2019,*) turbine%T_ave,turbine%P_ave,turbine%Torque_ave
      close(2019)

    end subroutine actuator_line_turbine_write_statistics

end module actuator_line_write_output
