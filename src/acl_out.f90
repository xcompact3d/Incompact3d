!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

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
      use param, only: itime, dt

      implicit none
      type(ActuatorLineType), intent(in) :: act_line
      character(len=100), intent(in) :: dir
      character(len=22) :: Format
      integer :: ielem
      logical :: exists

      do ielem=1,act_line%NElem
          inquire (File='loads/'//trim(act_line%name)//'_element_'//trim(int2str(ielem))//'.load', exist=exists)
          if (.not. exists) then
            call system('mkdir -p loads 2> /dev/null')
            open(2017,File='loads/'//trim(act_line%name)//'_element_'//trim(int2str(ielem))//'.load')
            write(2017,*) 'iteration,time,X,Y,Z,rdist/R,AOA,adot,RE,ur,CL,CD,CM25,Cn,Ct,Fn,Ft,F1'
            write(2017,*) '[-],[s],[m],[m],[m],[-],[deg],[deg/s],[-],[m/s],[-],[-],[-],[-],[-],[N],[N],[-]'
            close(2017)
          end if
          open(2017,File='loads/'//trim(act_line%name)//'_element_'//trim(int2str(ielem))//'.load', position="append", status="old", action="write")
          Format="(I5,A,20(E14.7,A))"
          write(2017,Format) itime,',',itime*dt,',',act_line%PEx(ielem),',',act_line%PEy(ielem),',',act_line%PEz(ielem),',',act_line%ERdist(ielem)/act_line%L,',',act_line%EAOA(ielem)*180/pi,',',act_line%EAOAdot(ielem)*pi/180,',',act_line%ERE(ielem),',',act_line%EUr(ielem),',',act_line%ECL(ielem),',',act_line%ECD(ielem),',',act_line%ECM(ielem),',',act_line%ECN(ielem),',',act_line%ECT(ielem),',',act_line%EFn(ielem),',',act_line%EFt(ielem),',',act_line%EEndeffects_factor(ielem)
          close(2017)
      enddo

    end subroutine actuator_line_element_write_output
    
    !*******************************************************************************
    !
    subroutine dynamic_stall_write_output(act_line,dir)
    !
    !*******************************************************************************
      use param, only: itime, dt

      implicit none
      type(ActuatorLineType), intent(in) :: act_line
      character(len=100), intent(in) :: dir
      character(len=22) :: Format
      integer :: ielem
      logical :: exists

      do ielem=1,act_line%NElem
          inquire (File='loads/'//trim(act_line%name)//'_element_'//trim(int2str(ielem))//'.dynstall', exist=exists)
          if (.not. exists) then
              call system('mkdir -p loads 2> /dev/null')
              open(2018,File='loads/'//trim(act_line%name)//'_element_'//trim(int2str(ielem))//'.dynstall')
              write(2018,*) 'iteration,time,rdist/R,pitch,AOA,f'
              close(2018)
          end if
          open(2018,File='loads/'//trim(act_line%name)//'_element_'//trim(int2str(ielem))//'.dynstall', position="append", status="old", action="write")
          Format="(I5,A,5(E14.7,A))"
          write(2018,Format)itime,',',itime*dt,',',act_line%ERdist(ielem)/act_line%L,',',act_line%Epitch(ielem)*180/pi,',',act_line%EAOA(ielem)*180/pi,',',act_line%EDynstall(ielem)%fprime
          close(2018)
      enddo

    end subroutine dynamic_stall_write_output
    
    !*******************************************************************************
    !
    subroutine actuator_line_turbine_write_output(turbine,dir)
    !
    !*******************************************************************************

      use param, only: itime, dt

      implicit none
      type(TurbineType), intent(in) :: turbine
      character(len=100), intent(in) :: dir
      character(len=22) :: Format
      logical :: exists

      !file that uses dir will be different each time step - dont use dir!
      inquire (File=trim(turbine%name)//'.perf', exist=exists)
      if (.not. exists) then
          open(2016,File=trim(turbine%name)//'.perf')
          write(2016,*) 'Iteration, Time, Number of Revs, GeneratorSpeed, GeneratorTorque, BladePitch1, BladePitch2, BladePitch3, Omega, DOmega, Ux, Uy, Uz, Thrust, Torque, Power'
          write(2016,*) '[-], [s], [-], [rad/s], [N m], [deg], [deg], [deg], [rad/s], [rad/s], [m/s], [m/s], [m/s], [N], [N m], [W]'
          close(2016)
      end if
      open(2016, File=trim(turbine%name)//'.perf', position="append", status="old", action="write")
      Format="(i5,A, 15(E14.7,A))"
      write(2016,Format) itime,',',itime*dt,',',turbine%AzimAngle/(2*pi),',',turbine%controller%GenSpeed,',',&
              turbine%controller%GenTrq,',',turbine%controller%PitCom(1)*180.0_mytype/pi,',',&
              turbine%controller%PitCom(2)*180.0_mytype/pi,',',turbine%controller%PitCom(3)*180.0_mytype/pi,',',&
              turbine%angularVel,',',turbine%deltaOmega,',',turbine%Ux_upstream,',',turbine%Uy_upstream,',',&
              turbine%Uz_upstream,',',turbine%Thrust,',',turbine%Torque,',',turbine%Power
      close(2016)

    end subroutine actuator_line_turbine_write_output   
    
    !*******************************************************************************
    !
    subroutine actuator_line_turbine_write_statistics(turbine,dir)
    !
    !*******************************************************************************

      use param, only: itime, dt

      implicit none
      type(TurbineType), intent(in) :: turbine
      character(len=100), intent(in) :: dir
      logical exists

      inquire (File=trim(turbine%name)//'.stat', exist=exists)
      if (.not. exists) then
          open(2019, File=trim(turbine%name)//'.stat')
          write(2019, *) 'Iteration Time T_ave, P_ave, Torque_ave'
          close(2019)
      end if
      open(2019, File=trim(turbine%name)//'.stat', position="append", status="old", action="write")
      write(2019,*) itime,itime*dt,turbine%T_ave,turbine%P_ave,turbine%Torque_ave
      close(2019)

    end subroutine actuator_line_turbine_write_statistics

end module actuator_line_write_output
