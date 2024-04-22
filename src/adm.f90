!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module actuator_disc_model

    use decomp_2d_mpi, only: nrank, nproc
    use decomp_2d_constants, only: mytype, real_type
    use actuator_line_model_utils 
    use airfoils

    implicit none
    
type ActuatorDiscType
    integer :: ID                       ! Actuator disc ID
    real(mytype) :: D                   ! Actuator disc diameter 
    real(mytype) :: COR(3)              ! Centre of rotation
    real(mytype) :: RotN(3)             ! Axis of rotation
    real(mytype) :: C_T                 ! Thrust coefficient
    real(mytype) :: alpha               ! Induction coefficient
    real(mytype) :: Udisc               ! Disc-averaged velocity
    real(mytype) :: UF                  ! Disc-averaged velocity -- filtered
    real(mytype) :: YawAng              ! Rotor yaw angle (in degrees, with respect to y-axis)
    real(mytype) :: TiltAng             ! Rotor tilt angle (in degrees, with respect to z-axis)
    real(mytype) :: Area
    real(mytype) :: Power
    real(mytype) :: Thrust
    real(mytype) :: Udisc_ave=0.0_mytype
    real(mytype) :: Power_ave=0.0_mytype
    real(mytype) :: Thrust_ave=0.0_mytype
end type ActuatorDiscType

    type(ActuatorDiscType), allocatable, save :: ActuatorDisc(:)
    integer, save :: Nad                ! Number of the actuator disc turbines 

contains

    !*******************************************************************************
    !
    subroutine actuator_disc_model_init(Ndiscs,admCoords)
    !
    !*******************************************************************************
        
      use var, only: GammaDisc
      use param, only: dx,dy,dz,istret,irestart,itime,iturboutput,initstat
      use decomp_2d
      use decomp_2d_io
      use MPI

      implicit none
      integer, intent(in) :: Ndiscs
      character(len=100), intent(in) :: admCoords
      character(1000) :: ReadLine
      real(mytype) :: GammaDisc_tot,GammaDisc_partial
      integer :: idisc,i,j,k,ierr,code
      real :: temp
      real(mytype) :: xmesh,ymesh,zmesh,deltax,deltay,deltaz,deltan,deltar,disc_thick,hgrid,projected_x,projected_y,projected_z

      ! ADM not yet set-up for stretched grids
      if (istret/=0) then
         write(*,*) 'Simulation stopped: run with different options'
         call MPI_ABORT(MPI_COMM_WORLD,code,ierr); stop
      endif

      ! Specify the actuator discs
      Nad=Ndiscs
      if (nrank==0) then
         write(*,*) '==========================================================='
         write(*,*) 'The actuator disc model is enabled'
         write(*,*) 'Number of Actuator discs : ', Nad
         write(*,*) '==========================================================='
      endif

      if (Nad>0) then 
         ! Read the disc data
         allocate(ActuatorDisc(Nad))
         open(15,file=admCoords)
         read(15,*)
         do idisc=1,Nad
            actuatordisc(idisc)%ID=idisc
            read(15,'(A)') ReadLine
            read(Readline,*) ActuatorDisc(idisc)%COR(1),ActuatorDisc(idisc)%COR(2),ActuatorDisc(idisc)%COR(3),&
                             ActuatorDisc(idisc)%YawAng,ActuatorDisc(idisc)%TiltAng,ActuatorDisc(idisc)%D,&
                             ActuatorDisc(idisc)%C_T,ActuatorDisc(idisc)%alpha
            ! Remember: RotN(1)=cos(yaw_angle)*cos(tilt_angle), RotN(2)=sin(tilt_angle) and RotN(3)=sin(yaw_angle)
            ActuatorDisc(idisc)%RotN(1)=cos(ActuatorDisc(idisc)%YawAng*conrad)*cos(ActuatorDisc(idisc)%TiltAng*conrad)  
            ActuatorDisc(idisc)%RotN(2)=sin(ActuatorDisc(idisc)%TiltAng*conrad)
            ActuatorDisc(idisc)%RotN(3)=sin(ActuatorDisc(idisc)%YawAng*conrad)  
            ActuatorDisc(idisc)%Area=pi*(ActuatorDisc(idisc)%D**2._mytype)/4._mytype
         enddo
         close(15)
          
         ! Compute Gamma
         GammaDisc=0.0
         do idisc=1,Nad
            GammaDisc_partial=0.0
            GammaDisc_tot=0.0
            ! Define the disc thickness 
            hgrid=sqrt((dx*actuatordisc(idisc)%RotN(1))**2+&
                       (dy*actuatordisc(idisc)%RotN(2))**2+&
                       (dz*(-actuatordisc(idisc)%RotN(3)))**2)
            disc_thick=max(actuatordisc(idisc)%D/8.0,hgrid*1.5)
            do k=1,xsize(3)
               zmesh=(xstart(3)+k-1-1)*dz
               do j=1,xsize(2)
                  ymesh=(xstart(2)+j-1-1)*dy
                  do i=1,xsize(1)
                     xmesh=(xstart(1)+i-1-1)*dx
                     ! Compute distance of mesh node to centre of disc 
                     ! [Quick fix to simulate a plate: set deltaz to zero and modify area]
                     deltax=xmesh-actuatordisc(idisc)%COR(1)
                     deltay=ymesh-actuatordisc(idisc)%COR(2)
                     deltaz=zmesh-actuatordisc(idisc)%COR(3)
                     
                     ! Projection of the vector (xmesh-COR(1), ymesh-COR(2), zmesh-COR(3)) onto the plane's normal vector (ROT(1), ROT(2), ROT(3))
                     ! deltan: distance from the point (xmesh,ymesh,zmesh) to the wind turbine plane along the plane's normal vector
                     deltan = deltax*actuatordisc(idisc)%RotN(1)+deltay*actuatordisc(idisc)%RotN(2)-deltaz*actuatordisc(idisc)%RotN(3)
                     
                     ! projected_: coordinates of the projected point
                     projected_x = xmesh-deltan*  actuatordisc(idisc)%RotN(1)
                     projected_y = ymesh-deltan*  actuatordisc(idisc)%RotN(2)
                     projected_z = zmesh-deltan*(-actuatordisc(idisc)%RotN(3))
                     
                     ! deltar: distance between the wind turbine's centre and the projected point
                     deltar = sqrt((projected_x-actuatordisc(idisc)%COR(1))**2+&
                                   (projected_y-actuatordisc(idisc)%COR(2))**2+&
                                   (projected_z-actuatordisc(idisc)%COR(3))**2)
                    
                     ! Compute Gamma [smearing using super-Gaussian functions, see also King et al. (2017) Wind Energ. Sci., 2, 115-131]
                     GammaDisc(i,j,k,idisc)=exp( -(deltan/(disc_thick/2.0))**2.0 -(deltar/(actuatordisc(idisc)%D/2.0))**8.0 )
                     GammaDisc_partial=GammaDisc_partial + GammaDisc(i,j,k,idisc) 
                  enddo
               enddo
            enddo
         
            call MPI_ALLREDUCE(GammaDisc_partial,GammaDisc_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
            GammaDisc(:,:,:,idisc)=GammaDisc(:,:,:,idisc)/GammaDisc_tot
            !if (nrank==0) then
            !   write(*,*) 'Disc thickness :', disc_thick
            !   write(*,*) 'Total Gamma volume: ', GammaDisc_tot*dx*dy*dz
            !   write(*,*) 'Approx. disc volume: ', actuatordisc(idisc)%Area*disc_thick
            !endif
         enddo

         if (irestart==1) then
             do idisc=1,Nad
                 open(15,File='disc'//trim(int2str(idisc))//'.adm', position="append", status='old', action='read')
                 backspace(15)
                 read(15,*) temp,temp,actuatordisc(idisc)%UF,actuatordisc(idisc)%Power,actuatordisc(idisc)%Thrust,actuatordisc(idisc)%Udisc_ave,actuatordisc(idisc)%Power_ave,actuatordisc(idisc)%Thrust_ave,temp,temp,temp
                 actuatordisc(idisc)%Udisc_ave=actuatordisc(idisc)%Udisc_ave*(itime-initstat+1)
                 actuatordisc(idisc)%Power_ave=actuatordisc(idisc)%Power_ave*(itime-initstat+1)
                 actuatordisc(idisc)%Thrust_ave=actuatordisc(idisc)%Thrust_ave*(itime-initstat+1)
                 close(15)
             end do
         endif

      endif

      return

    end subroutine actuator_disc_model_init

    !*******************************************************************************
    !
    subroutine actuator_disc_model_compute_source(ux1,uy1,uz1)
    !
    !*******************************************************************************
        
      use decomp_2d, only: xsize
      use MPI
      use param, only: dx, dy, dz, dt, itime, initstat, rho_air, T_relax, dBL, ustar
      use var, only: Fdiscx, Fdiscy, Fdiscz, GammaDisc
        
      implicit none
      real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1       
      real(mytype) :: uave, CTprime, alpha_relax
      real(mytype), allocatable, dimension(:) :: Udisc_partial
      integer :: i, j, k, idisc, ierr

      ! Compute disc-averaged velocity
      allocate(Udisc_partial(Nad))
     
      do idisc=1,Nad
         uave=0.
         do k=1,xsize(3)
            do j=1,xsize(2)
               do i=1,xsize(1)
                  ! Take the inner product to compute the rotor-normal velocity
                  uave=uave+GammaDisc(i,j,k,idisc)*(&  
                       ux1(i,j,k)*  actuatordisc(idisc)%RotN(1)+&
                       uy1(i,j,k)*  actuatordisc(idisc)%RotN(2)+&
                       uz1(i,j,k)*(-actuatordisc(idisc)%RotN(3)))
               enddo
            enddo
         enddo
         Udisc_partial(idisc)=uave
      enddo
      call MPI_ALLREDUCE(Udisc_partial,actuatordisc%Udisc,Nad,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)

      deallocate(Udisc_partial)

      ! Time relaxation -- low pass filter
      if (itime==1.or.T_relax<0) then
         do idisc=1,Nad
            actuatordisc(idisc)%UF=actuatordisc(idisc)%Udisc
         enddo
      else
         alpha_relax=(dt/T_relax)/(1.+dt/T_relax)
         do idisc=1,Nad
            actuatordisc(idisc)%UF=alpha_relax*actuatordisc(idisc)%Udisc+(1.-alpha_relax)*actuatordisc(idisc)%UF
         enddo
      endif
        
      ! Compute the forces
      Fdiscx=0.
      Fdiscy=0.
      Fdiscz=0.
      do idisc=1,Nad
         ! Compute local thrust coefficient
         CTprime=actuatordisc(idisc)%C_T/(1-actuatordisc(idisc)%alpha)**2.
         ! Compute power and thrust (a power coefficient may be considered here)
         actuatordisc(idisc)%Thrust=0.5_mytype*rho_air*CTprime*actuatordisc(idisc)%UF**2.0_mytype*actuatordisc(idisc)%Area
         actuatordisc(idisc)%Power =actuatordisc(idisc)%Thrust*actuatordisc(idisc)%UF
         ! Distribute thrust force
         Fdiscx(:,:,:)=Fdiscx(:,:,:)-actuatordisc(idisc)%Thrust*GammaDisc(:,:,:,idisc)*actuatordisc(idisc)%RotN(1)/(dx*dy*dz)
         Fdiscy(:,:,:)=Fdiscy(:,:,:)-actuatordisc(idisc)%Thrust*GammaDisc(:,:,:,idisc)*actuatordisc(idisc)%RotN(2)/(dx*dy*dz)
         Fdiscz(:,:,:)=Fdiscz(:,:,:)-actuatordisc(idisc)%Thrust*GammaDisc(:,:,:,idisc)*(-actuatordisc(idisc)%RotN(3))/(dx*dy*dz)
      enddo
        
      ! Compute statistics
      if (itime>=initstat) then
         do idisc=1,Nad
            actuatordisc(idisc)%Udisc_ave=actuatordisc(idisc)%Udisc_ave+actuatordisc(idisc)%UF
            actuatordisc(idisc)%Power_ave=actuatordisc(idisc)%Power_ave+actuatordisc(idisc)%Power
            actuatordisc(idisc)%Thrust_ave=actuatordisc(idisc)%Thrust_ave+actuatordisc(idisc)%Thrust 
         enddo
      endif

      return

    end subroutine actuator_disc_model_compute_source
 
    !*******************************************************************************
    !
    subroutine actuator_disc_model_write_output(dump_no)
    !
    !*******************************************************************************

      use param, only: itime, initstat, dt

      implicit none
      integer, intent(in) :: dump_no
      integer :: idisc

      if (Nad>0) then
         do idisc=1, Nad
             if (dump_no==1) then
                 open(2020,File='disc'//trim(int2str(idisc))//'.adm')
                 write(2020,*) 'i, Time, UF, Power, Thrust, Udisc_ave, Power_ave, Thrust_ave, alpha, YawAng, TiltAng'
                 close(2020)
             endif
             open(2020,File='disc'//trim(int2str(idisc))//'.adm', position="append", status="old", action="write")
             write(2020,*)itime,&
                          itime*dt,&
                          actuatordisc(idisc)%UF,&
                          actuatordisc(idisc)%Power,&
                          actuatordisc(idisc)%Thrust,&
                          actuatordisc(idisc)%Udisc_ave/(itime-initstat+1),&
                          actuatordisc(idisc)%Power_ave/(itime-initstat+1),&
                          actuatordisc(idisc)%Thrust_ave/(itime-initstat+1),&
                          actuatordisc(idisc)%alpha,&
                          actuatordisc(idisc)%YawAng,&
                          actuatordisc(idisc)%TiltAng
             close(2020)
         enddo
      endif

      return 

    end subroutine actuator_disc_model_write_output

end module actuator_disc_model
