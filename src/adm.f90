module actuator_disc_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank
    use actuator_line_model_utils 
    use airfoils

    implicit none
    
    type ActuatorDiscType
        integer :: ID                       ! Actuator disk ID
        real(mytype) :: D                   ! Actuator disk diameter 
        real(mytype) :: COR(3)              ! Center of Rotation
        real(mytype) :: RotN(3)             ! axis of rotation
        real(mytype) :: C_T                 ! Thrust coefficient
        real(mytype) :: alpha               ! Induction coefficient
        real(mytype) :: Udisc               ! Disc-averaged velocity
        real(mytype) :: UF                  ! Disc-averaged velocity -- filtered
        real(mytype) :: UF_prev             ! Disc-averaged velocity -- filtered previous
        real(mytype) :: Power
        real(mytype) :: Thrust
        real(mytype) :: Udisc_ave=0.0_mytype
        real(mytype) :: Power_ave=0.0_mytype
        real(mytype) :: Thrust_ave=0.0_mytype
    end type ActuatorDiscType

    type(ActuatorDiscType), allocatable, save :: ActuatorDisc(:)
    integer,save :: Nad ! Number of the actuator disk turbines 

contains

    subroutine actuator_disc_model_init(Ndiscs,admCoords,C_T,aind)
        
        USE var, only: Fdiscx, Fdiscy, Fdiscz, GammaDisc, yp
        use param, only: dx,dy,dz,istret
        USE decomp_2d
        USE decomp_2d_io
        USE MPI

        implicit none
        integer,intent(in) :: Ndiscs
        character(len=100),intent(in) :: admCoords
        character(1000) :: ReadLine
        real(mytype) :: C_T, aind
        integer :: idisc,i,j,k
        real(mytype) :: xmesh,ymesh,zmesh,deltax,deltay,deltaz,deltar,dr

        !### Specify the actuator discs
        Nad=Ndiscs
        if (nrank==0) then
          print *, "==========================================================="
          print *, "The actuator disc model is enabled"
          print *, 'Number of Actuator discs : ', Nad
        endif

        if (Nad>0) then 
          allocate(ActuatorDisc(Nad))
          open(15,file=admCoords)
          do idisc=1,Nad
            actuatordisc(idisc)%ID=idisc
            read(15,'(A)') ReadLine
            read(Readline,*) ActuatorDisc(idisc)%COR(1),ActuatorDisc(idisc)%COR(2),ActuatorDisc(idisc)%COR(3),ActuatorDisc(idisc)%RotN(1),ActuatorDisc(idisc)%RotN(2),ActuatorDisc(idisc)%RotN(3),ActuatorDisc(idisc)%D 
            ActuatorDisc(idisc)%C_T=C_T
            ActuatorDisc(idisc)%alpha=aind
          enddo
          close(15)
          
          ! Compute Gamma
          GammaDisc=0.
          do idisc=1,Nad
            do k=1,xsize(3)
              zmesh=(xstart(3)+k-1-1)*dz
              do j=1,xsize(2)
                if (istret.eq.0) ymesh=(xstart(2)+j-1-1)*dy
                if (istret.ne.0) ymesh=yp(xstart(2)+j)
                do i=1,xsize(1)
                  xmesh=(xstart(1)+i-1-1)*dx
                  deltax=abs(xmesh-actuatordisc(idisc)%COR(1))
                  deltay=abs(ymesh-actuatordisc(idisc)%COR(2))
                  deltaz=abs(zmesh-actuatordisc(idisc)%COR(3))
                  deltar=sqrt(deltay**2.+deltaz**2.)
                  dr=sqrt(dy**2.+dz**2.)

                  if(deltax>dx/2.) then
                    GammaDisc(i,j,k,idisc)=0.
                  elseif (deltax<=dx/2.) then
                    if(deltar<=actuatordisc(idisc)%D/2.0) then
                      GammaDisc(i,j,k,idisc)=1.
                    elseif(deltar>actuatordisc(idisc)%D/2..and. deltar<=actuatordisc(idisc)%D/2.+dr) then
                      GammaDisc(i,j,k,idisc)=1.-(deltar-actuatordisc(idisc)%D/2.)/dr
                    else
                      GammaDisc(i,j,k,idisc)=0.
                    endif
                  endif

                enddo
              enddo
            enddo
          enddo
        endif

        return

    end subroutine actuator_disc_model_init

    subroutine actuator_disc_model_compute_source(ux1,uy1,uz1)
        
        use decomp_2d, only: mytype, nproc, xstart, xend, xsize, update_halo, decomp_2d_abort
        use MPI
        use param, only: dx,dy,dz,istret,dt,itime,initstat,rho_air,dBL,ustar,ifirst
        use var, only: FDiscx, FDiscy, FDiscz, GammaDisc, yp
        
        implicit none
        real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1       
        real(mytype) :: xmesh, ymesh,zmesh,deltax,deltay,deltaz,deltar,DiscsTotalArea
        real(mytype) :: uave,CTprime, T_relax, alpha_relax, Uinf, CTave,Ratio,sumforce,Sumforce_partial
        real(mytype), allocatable, dimension(:) :: Udisc_partial
        integer,allocatable, dimension(:) :: counter, counter_total
        integer :: i,j,k, idisc, code

        ! Compute disc-averaged velocity
        allocate(Udisc_partial(Nad))
        allocate(counter(Nad))
        allocate(counter_total(Nad))
     
        do idisc=1,Nad
          uave=0.
          counter(idisc)=0
          counter_total(idisc)=0.
          do k=1,xsize(3)
            zmesh=(xstart(3)+k-1-1)*dz 
            do j=1,xsize(2)
              if (istret.eq.0) ymesh=(xstart(2)+j-1-1)*dy
              if (istret.ne.0) ymesh=yp(xstart(2)+j-1)
              do i=1,xsize(1)
                xmesh=(i-1)*dx
                deltax=abs(xmesh-actuatordisc(idisc)%COR(1))
                deltay=abs(ymesh-actuatordisc(idisc)%COR(2))
                deltaz=abs(zmesh-actuatordisc(idisc)%COR(3))
                deltar=sqrt(deltay**2.+deltaz**2.)
                if(deltar<=actuatordisc(idisc)%D/2. .and. deltax<=dx/2.) then
                  ! Take the inner product to compute the rotor-normal velocity
                  uave=uave+&!sqrt(ux1(i,j,k)**2.+uy1(i,j,k)**2.+uz1(i,j,k)**2.)
                       ux1(i,j,k)*actuatordisc(idisc)%RotN(1)+&
                       uy1(i,j,k)*actuatordisc(idisc)%RotN(2)+&
                       uz1(i,j,k)*actuatordisc(idisc)%RotN(3)
                  counter(idisc)=counter(idisc)+1
                endif
              enddo
            enddo
          enddo
          Udisc_partial(idisc)=uave
        enddo

        call MPI_ALLREDUCE(Udisc_partial,actuatordisc%Udisc,Nad,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")
        call MPI_ALLREDUCE(counter,counter_total,Nad,MPI_INTEGER,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

        do idisc=1,Nad
          if (counter_total(idisc)==0) then
            print *, 'counter=0 for the disc --- something is wrong with the disc number: ', idisc
            stop
          endif
          actuatordisc(idisc)%Udisc=actuatordisc(idisc)%Udisc/counter_total(idisc)
        enddo
        
        deallocate(Udisc_partial,counter,counter_total)

        ! Time relaxation -- low pass filter
        T_relax= 0.27*dBL/ustar !5.0, 5.0*dt
        alpha_relax=(dt/T_relax)/(1.+dt/T_relax)
        if (itime==ifirst) then
          do idisc=1,Nad
            actuatordisc(idisc)%UF=actuatordisc(idisc)%Udisc
            !actuatordisc(idisc)%UF_prev=actuatordisc(idisc)%Udisc
          enddo
        else
          do idisc=1,Nad
            actuatordisc(idisc)%UF=alpha_relax*actuatordisc(idisc)%Udisc+(1.-alpha_relax)*actuatordisc(idisc)%UF
            !actuatordisc(idisc)%UF=alpha_relax*actuatordisc(idisc)%UF+(1.-alpha_relax)*actuatordisc(idisc)%UF_prev
            !actuatordisc(idisc)%UF_prev=actuatordisc(idisc)%Udisc
          enddo
        endif
        
        ! Compute the forces
        Fdiscx=0.
        do idisc=1,Nad
          CTprime=actuatordisc(idisc)%C_T/(1-actuatordisc(idisc)%alpha)**2.
          Fdiscx(:,:,:)=Fdiscx(:,:,:)-0.5_mytype*rho_air*CTprime*actuatordisc(idisc)%UF**2.0_mytype*GammaDisc(:,:,:,idisc)/dx
          ! Compute power and thrust
          actuatordisc(idisc)%Power =0.5_mytype*rho_air*CTprime*actuatordisc(idisc)%UF**2.0_mytype*actuatordisc(idisc)%UF&
                                    *pi*actuatordisc(idisc)%D**2._mytype/4._mytype
          actuatordisc(idisc)%Thrust=0.5_mytype*rho_air*CTprime*actuatordisc(idisc)%UF**2.0_mytype&
                                    *pi*actuatordisc(idisc)%D**2._mytype/4._mytype
        enddo
        Fdiscy(:,:,:)=0.
        Fdiscz(:,:,:)=0.
        
        !if(iverifyadm.eq.1) then
        !  sumforce=0.
        !  do k=1,xsize(3)
        !  do j=1,xsize(2)
        !  do i=1,xsize(1)
        !    sumforce_partial=sumforce_partial+Fdiscx(i,j,k)
        !  enddo
        !  enddo
        !  enddo
        !  call MPI_ALLREDUCE(Sumforce_partial,sumforce,1,MPI_REAL8,MPI_SUM, &
        !       MPI_COMM_WORLD,ierr)
        !  do idisc=1,Nad
        !    DiscsTotalArea=DiscsTotalArea+pi/4.0_mytype*actuatordisc(idisc)%D**2.0_mytype
        !  enddo
        !  CTave=sum(actuatordisc%C_T)/Nad
        !  Uinf=sum(actuatordisc%UF)/Nad
        !  Ratio=sumforce*dx*dy*dz/(-0.5*CTprime*Uinf**2*DiscsTotalArea)
        !    if (nrank==0) print *, 'ADM verification ratio', Ratio
        !    if (nrank==0) print *, 'Thrust', sumforce*dx*dy*dz
        !endif
        
        ! Compute the average values after spin_up time
        if(itime>=initstat) then
          do idisc=1,Nad
            actuatordisc(idisc)%Udisc_ave=actuatordisc(idisc)%Udisc_ave+actuatordisc(idisc)%UF
            actuatordisc(idisc)%Power_ave=actuatordisc(idisc)%Power_ave+actuatordisc(idisc)%Power 
            actuatordisc(idisc)%Thrust_ave=actuatordisc(idisc)%Thrust_ave+actuatordisc(idisc)%Thrust 
          enddo
        endif

        return
    end subroutine actuator_disc_model_compute_source

 
    subroutine actuator_disc_model_write_output(dump_no)

        implicit none
        integer,intent(in) :: dump_no
        integer :: idisc
        character(len=300) :: dir, Format

        if (Nad>0) then
          open(2020,File='discs_time'//trim(int2str(dump_no))//'.adm')
          write(2020,*) 'Udisc, CT, Power, Thrust, Udisc_ave, Power_ave, Thrust_ave'
          Format="(7(E14.7,A))"
          do idisc=1,Nad
            write(2020,Format) actuatordisc(idisc)%UF,',',actuatordisc(idisc)%C_T,',',actuatordisc(idisc)%Power,',',actuatordisc(idisc)%Thrust,',',actuatordisc(idisc)%Udisc_ave,',',actuatordisc(idisc)%Power_ave,',',actuatordisc(idisc)%Thrust_ave
          end do
          close(2020)
        endif

        return 
    end subroutine actuator_disc_model_write_output
end module actuator_disc_model
