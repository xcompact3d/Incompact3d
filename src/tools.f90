!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module tools

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  use utilities

  implicit none

  logical, save :: adios2_restart_initialised = .false.

  character(len=*), parameter :: io_restart = "restart-io"
  character(len=*), parameter :: resfile = "checkpoint"
  character(len=*), parameter :: resfile_old = "checkpoint.old"
  character(len=*), parameter :: io_ioflow = "in-outflow-io"
  
  private

  public :: test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, &
       apply_spatial_filter, read_inflow, append_outflow, write_outflow, init_inflow_outflow, &
       compute_cfldiff, compute_cfl, &
       rescale_pressure, mean_plane_x, mean_plane_y, mean_plane_z

contains
  !##################################################################
  !##################################################################
  subroutine test_scalar_min_max(phi)

    use variables
    use param
    use var
    use mpi

    implicit none

    integer :: code,ierror,i,j,k,is,jglob
    real(mytype) :: phimax,phimin,phimax1,phimin1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(2,numscalar) :: phimaxin,phimaxout

    do is=1, numscalar

      ta1(:,:,:) = phi(:,:,:,is)
      ! ibm
      if (iibm > 0) then
        ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
      endif

      phimax=-1609._mytype
      phimin=1609._mytype
      phimax = maxval(ta1(:,:,:))
      phimin =-minval(ta1(:,:,:))
      phimaxin(:,is) =  (/phimin, phimax /)
    enddo

    call MPI_REDUCE(phimaxin,phimaxout,numscalar*2,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    do is=1,numscalar
      if (nrank == 0) then
        phimin1 = -phimaxout(1,is)
        phimax1 =  phimaxout(2,is)

        write(*,*) 'Phi'//char(48+is)//' min max=', real(phimin1,4), real(phimax1,4)

        if (abs(phimax1) > 100._mytype) then !if phi control turned off
           write(*,*) 'Scalar diverged! SIMULATION IS STOPPED!'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
        endif
      endif
    enddo

    return
  end subroutine test_scalar_min_max
  !##################################################################
  !##################################################################
  subroutine test_speed_min_max(ux,uy,uz)

    use variables
    use param
    use var
    use mpi

    implicit none

    integer :: code,ierror,i,j,k
    real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin
    real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(6) :: umaxin, umaxout

    if (iibm > 0) then
       ux(:,:,:) = (one - ep1(:,:,:)) * ux(:,:,:)
       uy(:,:,:) = (one - ep1(:,:,:)) * uy(:,:,:)
       uz(:,:,:) = (one - ep1(:,:,:)) * uz(:,:,:)
    endif

    uxmax=-1609.;uymax=-1609.;uzmax=-1609.;uxmin=1609.;uymin=1609.;uzmin=1609.
    !
    ! More efficient version
    uxmax=maxval(ux)
    uymax=maxval(uy)
    uzmax=maxval(uz)
    uxmin=-minval(ux)
    uymin=-minval(uy)
    uzmin=-minval(uz)

    umaxin = (/uxmax, uymax, uzmax, uxmin, uymin, uzmin/)
    call MPI_REDUCE(umaxin,umaxout,6,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    uxmax1= umaxout(1)
    uymax1= umaxout(2)
    uzmax1= umaxout(3)
    uxmin1=-umaxout(4)
    uymin1=-umaxout(5)
    uzmin1=-umaxout(6)

    if (nrank == 0) then

       write(*,*) 'U,V,W min=',real(uxmin1,4),real(uymin1,4),real(uzmin1,4)
       write(*,*) 'U,V,W max=',real(uxmax1,4),real(uymax1,4),real(uzmax1,4)
       !print *,'CFL=',real(abs(max(uxmax1,uymax1,uzmax1)*dt)/min(dx,dy,dz),4)

       if((abs(uxmax1)>=onehundred).or.(abs(uymax1)>=onehundred).OR.(abs(uzmax1)>=onehundred)) then
         write(*,*) 'Velocity diverged! SIMULATION IS STOPPED!'
         call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
         stop
       endif

    endif

    return
  end subroutine test_speed_min_max
  !##################################################################
  !##################################################################
  subroutine simu_stats(iwhen)

    use simulation_stats
    use var
    use MPI

    implicit none

    integer :: iwhen

    if (iwhen == 1) then !AT THE START OF THE SIMULATION
       tstart=zero
       time1=zero
       trank=zero
       tranksum=zero
       ttotal=zero
       call cpu_time(tstart)
    else if (iwhen == 2) then !AT THE START OF A TIME STEP
       if (nrank == 0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
          call cpu_time(time1)
          write(*,*) '==========================================================='
          write(*,"(' Time step =',i7,'/',i7,', Time unit =',F12.4)") itime,ilast,t
       endif
    else if ((iwhen == 3).and.(itime > ifirst)) then !AT THE END OF A TIME STEP
       if (nrank == 0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
          call cpu_time(trank)
          if (nrank==0) write(*,*) 'Time for this time step (s):',real(trank-time1)
          telapsed = (trank-tstart)/threethousandsixhundred
          tremaining  = telapsed*(ilast-itime)/(itime-ifirst)
          write(*,"(' Remaining time:',I8,' h ',I2,' min')") int(tremaining), int((tremaining-int(tremaining))*sixty)
          write(*,"(' Elapsed time:  ',I8,' h ',I2,' min')") int(telapsed), int((telapsed-int(telapsed))*sixty)
       endif
    else if (iwhen == 4) then !AT THE END OF THE SIMULATION
       call cpu_time(trank)
       ttotal=trank-tstart
       if (nrank == 0) then
          write(*,*) '==========================================================='
          write(*,*) '                                                           '
          write(*,*) 'Good job! Xcompact3d finished successfully!                '
          write(*,*) '                                                           '
          write(*,*) '2DECOMP with p_row*p_col=',p_row,p_col
          write(*,*) '                                                           '
          write(*,*) 'nx*ny*nz=',nx*ny*nz
          write(*,*) 'nx,ny,nz=',nx,ny,nz
          write(*,*) 'dx,dy,dz=',dx,dy,dz
          write(*,*) '                                                           '
          write(*,*) 'Averaged time per step (s):',real(ttotal/(ilast-(ifirst-1)),4)
          write(*,*) 'Total wallclock (s):',real(ttotal,4)
          write(*,*) 'Total wallclock (m):',real(ttotal/sixty,4)
          write(*,*) 'Total wallclock (h):',real(ttotal/threethousandsixhundred,4)
          write(*,*) '                                                           '
       endif
    endif

  end subroutine simu_stats
  !##############################################################################
    !!
    !!  SUBROUTINE: restart
    !! DESCRIPTION: reads or writes restart file
    !!
    !!      AUTHOR: ?
    !!    MODIFIED: Kay Schäfer
    !!
  !##############################################################################
  subroutine restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,iresflg)

    use decomp_2d_io
    use variables
    use param
    use MPI
    use mod_stret, only : beta
    use navier, only : gradp
    use mhd, only : mhd_equation,Bm,dBm
    use particle, only : particle_checkpoint

    implicit none

    integer :: i,j,k,iresflg,nzmsize,ierror,is,it,code
    integer :: ierror_o=0 !error to open sauve file during restart
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: mu1
    real(mytype) :: xdt,tfield,y
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods
    logical :: fexists
    character(len=30) :: filename, filestart
    character(len=32) :: fmt2,fmt3,fmt4
    character(len=7) :: fmt1
    character(len=80) :: varname
    NAMELIST /Time/ tfield, itime
    NAMELIST /NumParam/ nx, ny, nz, istret, beta, dt, itimescheme

    logical, save :: first_restart = .true.
    
    write(filename,"('restart',I7.7)") itime
    write(filestart,"('restart',I7.7)") ifirst-1

    if (iresflg == 1) then !Writing restart
       if (mod(itime, icheckpoint) /= 0) then
          return
       endif

       if (nrank==0) then
          write(*,*) '===========================================================<<<<<'
          write(*,*) 'Writing restart point ',filename !itime/icheckpoint
       endif

       if (adios2_restart_initialised) then
          ! First, rename old checkpoint in case of error
          if (validation_restart) call rename(resfile, resfile_old)
       end if
    end if

    if (.not. adios2_restart_initialised) then
       call init_restart_adios2()
       adios2_restart_initialised = .true.
    end if

    if (iresflg==1) then !write
       call decomp_2d_open_io(io_restart, resfile, decomp_2d_write_mode)
       call decomp_2d_start_io(io_restart, resfile)

       call decomp_2d_write_one(1,ux1,resfile,"ux",0,io_restart,reduce_prec=.false.)
       call decomp_2d_write_one(1,uy1,resfile,"uy",0,io_restart,reduce_prec=.false.)
       call decomp_2d_write_one(1,uz1,resfile,"uz",0,io_restart,reduce_prec=.false.)
       ! write previous time-step if necessary for AB2 or AB3
       if ((itimescheme==2).or.(itimescheme==3)) then
          call decomp_2d_write_one(1,dux1(:,:,:,2),resfile,"dux-2",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duy1(:,:,:,2),resfile,"duy-2",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duz1(:,:,:,2),resfile,"duz-2",0,io_restart,reduce_prec=.false.)
       end if
       ! for AB3 one more previous time-step
       if (itimescheme==3) then
          call decomp_2d_write_one(1,dux1(:,:,:,3),resfile,"dux-3",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duy1(:,:,:,3),resfile,"duy-3",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duz1(:,:,:,3),resfile,"duz-3",0,io_restart,reduce_prec=.false.)
       end if
       !
       call decomp_2d_write_one(3,pp3,resfile,"pp",0,io_restart,opt_decomp=phG,reduce_prec=.false.)
       !
       if (iscalar==1) then
          do is=1, numscalar
             write(varname,"('phi-',I2.2)") is
             call decomp_2d_write_one(1,phi1(:,:,:,is),resfile,trim(varname),0,io_restart,reduce_prec=.false.)
             ! previous time-steps
             if ((itimescheme==2).or.(itimescheme==3)) then ! AB2 or AB3
                write(varname,"('dphi-',I2.2,'-2')") is
                call decomp_2d_write_one(1,dphi1(:,:,:,2,is),resfile,trim(varname),0,io_restart,reduce_prec=.false.)
             end if
             !
             if (itimescheme==3) then ! AB3
               write(varname,"('dphi-',I2.2,'-3')") is
               call decomp_2d_write_one(1,dphi1(:,:,:,3,is),resfile,trim(varname),0,io_restart,reduce_prec=.false.)
             end if
          end do
       endif
       if (ilmn) then
          do is = 1, nrhotime
             write(varname,"('rho-',I2.2)") is
             call decomp_2d_write_one(1,rho1(:,:,:,is),resfile,trim(varname),0,io_restart,reduce_prec=.false.)
          enddo
          do is = 1, ntime
             write(varname,"('drho-',I2.2)") is
             call decomp_2d_write_one(1,drho1(:,:,:,is),resfile,trim(varname),0,io_restart,reduce_prec=.false.)
          enddo
          call decomp_2d_write_one(1,mu1(:,:,:),resfile,"mu",0,io_restart,reduce_prec=.false.)
       endif

       if (mhd_active .and. mhd_equation == 'induction') then
          call decomp_2d_write_one(1,Bm(:,:,:,1),resfile,'bx',0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,Bm(:,:,:,2),resfile,'by',0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,Bm(:,:,:,3),resfile,'bz',0,io_restart,reduce_prec=.false.)
          if (itimescheme==3) then
             call decomp_2d_write_one(1,dBm(:,:,:,1,2),resfile,'dbx-2',0,io_restart,reduce_prec=.false.)
             call decomp_2d_write_one(1,dBm(:,:,:,2,2),resfile,'dby-2',0,io_restart,reduce_prec=.false.)
             call decomp_2d_write_one(1,dBm(:,:,:,3,2),resfile,'dbz-2',0,io_restart,reduce_prec=.false.)
             call decomp_2d_write_one(1,dBm(:,:,:,1,3),resfile,'dbx-3',0,io_restart,reduce_prec=.false.)
             call decomp_2d_write_one(1,dBm(:,:,:,2,3),resfile,'dby-3',0,io_restart,reduce_prec=.false.)
             call decomp_2d_write_one(1,dBm(:,:,:,3,3),resfile,'dbz-3',0,io_restart,reduce_prec=.false.)
          endif
       endif

       call decomp_2d_end_io(io_restart, resfile)
       call decomp_2d_close_io(io_restart, resfile)

       ! Validate restart file then remove old file and update restart.info
       if (validation_restart) then
          if (validate_restart(resfile_old, resfile)) then
             call delete_filedir(resfile_old)
          else
             call decomp_2d_abort(1, "Writing restart - validation failed!")
          endif
       endif

       ! Write info file for restart - Kay Schäfer
       if (nrank == 0) then
          write(filename,"('restart.info')")
          write(fmt2,'("(A,I16)")')
          write(fmt3,'("(A,F16.4)")')
          write(fmt4,'("(A,F16.12)")')
          !
          open (111,file=filename,action='write',status='replace')
          write(111,'(A)')'!========================='
          write(111,'(A)')'&Time'
          write(111,'(A)')'!========================='
          write(111,fmt3) 'tfield=   ',t
          write(111,fmt2) 'itime=    ',itime
          write(111,'(A)')'/End'
          write(111,'(A)')'!========================='
          write(111,'(A)')'&NumParam'
          write(111,'(A)')'!========================='
          write(111,fmt2) 'nx=       ',nx
          write(111,fmt2) 'ny=       ',ny
          write(111,fmt2) 'nz=       ',nz
          write(111,fmt3) 'Lx=       ',xlx
          write(111,fmt3) 'Ly=       ',yly
          write(111,fmt3) 'Lz=       ',zlz
          write(111,fmt2) 'istret=   ',istret
          write(111,fmt4) 'beta=     ',beta
          write(111,fmt2) 'iscalar=  ',iscalar
          write(111,fmt2) 'numscalar=',numscalar
          write(111,'(A,I14)') 'itimescheme=',itimescheme
          write(111,fmt2) 'iimplicit=',iimplicit
          write(111,'(A)')'/End'
          write(111,'(A)')'!========================='
          close(111)
       end if

       if(particle_active) then
          call particle_checkpoint(mode='write',filename='checkpoint-particles')
       endif

    else
       if (nrank==0) then
         write(*,*)'==========================================================='
         write(*,*)'RESTART from file:', filestart
         write(*,*)'==========================================================='
       end if
       call decomp_2d_open_io(io_restart, resfile, decomp_2d_read_mode)
       call decomp_2d_start_io(io_restart, resfile)

       call decomp_2d_read_one(1,ux1,resfile,"ux",io_restart,reduce_prec=.false.)
       call decomp_2d_read_one(1,uy1,resfile,"uy",io_restart,reduce_prec=.false.)
       call decomp_2d_read_one(1,uz1,resfile,"uz",io_restart,reduce_prec=.false.)
       ! read previous time-step if necessary for AB2 or AB3
       if ((itimescheme==2).or.(itimescheme==3)) then ! AB2 or AB3
          call decomp_2d_read_one(1,dux1(:,:,:,2),resfile,"dux-2",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duy1(:,:,:,2),resfile,"duy-2",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duz1(:,:,:,2),resfile,"duz-2",io_restart,reduce_prec=.false.)
       end if
       ! for AB3 one more previous time-step
       if (itimescheme==3) then ! AB3
          call decomp_2d_read_one(1,dux1(:,:,:,3),resfile,"dux-3",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duy1(:,:,:,3),resfile,"duy-3",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duz1(:,:,:,3),resfile,"duz-3",io_restart,reduce_prec=.false.)
       end if
       !
       call decomp_2d_read_one(3,pp3,resfile,"pp",io_restart,opt_decomp=phG,reduce_prec=.false.)
       !
       if (iscalar==1) then
         do is=1, numscalar
            write(varname,"('phi-',I2.2)") is
            call decomp_2d_read_one(1,phi1(:,:,:,is),resfile,trim(varname),io_restart,reduce_prec=.false.)
           ! previous time-steps
           if ((itimescheme==2).or.(itimescheme==3)) then ! AB2 or AB3
             write(varname,"('dphi-',I2.2,'-2')") is
             call decomp_2d_read_one(1,dphi1(:,:,:,2,is),resfile,trim(varname),io_restart,reduce_prec=.false.)
           end if
           !
           if (itimescheme==3) then ! AB3
              write(varname,"('dphi-',I2.2,'-3')") is
              call decomp_2d_read_one(1,dphi1(:,:,:,3,is),resfile,trim(varname),io_restart,reduce_prec=.false.)
           end if
           ! ABL 
           if (itype==itype_abl) then
             do j=1,xsize(2)
               if (istret==0) y = (j + xstart(2)-1-1)*dy
               if (istret.ne.0) y = yp(j+xstart(2)-1)
               if (ibuoyancy==1) then
                 Tstat(j,1) = T_wall - (T_wall-T_top)*y/yly
               else
                 Tstat(j,1) = zero
               endif
             enddo
           endif
         end do
       endif
       if (ilmn) then
          do is = 1, nrhotime
             write(varname,"('rho-',I2.2)") is
             call decomp_2d_read_one(1,rho1(:,:,:,is),resfile,trim(varname),io_restart,reduce_prec=.false.)
          enddo
          do is = 1, ntime
             write(varname,"('drho-',I2.2)") is
             call decomp_2d_read_one(1,drho1(:,:,:,is),resfile,trim(varname),io_restart,reduce_prec=.false.)
          enddo
          call decomp_2d_read_one(1,mu1,resfile,"mu",io_restart,reduce_prec=.false.)
       end if

       if(mhd_active .and. mhd_equation == 'induction') then
          call decomp_2d_read_one(1,Bm(:,:,:,1),resfile,'bx',io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,Bm(:,:,:,2),resfile,'by',io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,Bm(:,:,:,3),resfile,'bz',io_restart,reduce_prec=.false.)
          if (itimescheme==3) then
             call decomp_2d_read_one(1,dBm(:,:,:,1,2),resfile,'dbx-2',io_restart,reduce_prec=.false.)
             call decomp_2d_read_one(1,dBm(:,:,:,2,2),resfile,'dby-2',io_restart,reduce_prec=.false.)
             call decomp_2d_read_one(1,dBm(:,:,:,3,2),resfile,'dbz-2',io_restart,reduce_prec=.false.)
             call decomp_2d_read_one(1,dBm(:,:,:,1,3),resfile,'dbx-3',io_restart,reduce_prec=.false.)
             call decomp_2d_read_one(1,dBm(:,:,:,2,3),resfile,'dby-3',io_restart,reduce_prec=.false.)
             call decomp_2d_read_one(1,dBm(:,:,:,3,3),resfile,'dbz-3',io_restart,reduce_prec=.false.)
          endif
       endif

       call decomp_2d_end_io(io_restart, resfile)
       call decomp_2d_close_io(io_restart, resfile)

       !! Read time of restart file
       write(filename,"(A)") 'restart.info'
       inquire(file=filename, exist=fexists)
       if (nrank==0) write(*,*) filename
       ! file exists???
       if (fexists) then
         open(111, file=filename)
         read(111, nml=Time)
         close(111)

         t0 = tfield
         itime0 = 0

       else
         t0 = zero
         itime0 = 0
       end if

       if(particle_active) call particle_checkpoint(mode='read')
       
    endif

    if (nrank==0) then
       if (ierror_o /= 0) then !Included by Felipe Schuch
          write(*,*) '==========================================================='
          write(*,*) 'Error: Impossible to read '//trim(filestart)
          write(*,*) '==========================================================='
          call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
       endif
    endif

    ! reconstruction of the dp/dx, dp/dy and dp/dz from pp3
    if (iresflg==0) then
       if (itimescheme <= 4) itr=1
       if (itimescheme == 5) itr=3
       if (itimescheme == 6) itr=5
       call gradp(px1,py1,pz1,pp3)
       if (nrank == 0) write(*,*) 'reconstruction pressure gradients done!'
    end if

    if (iresflg==1) then !Writing restart
       if (nrank==0) then
          write(fmt1,"(I7.7)") itime
          write(*,*) 'Restart point restart',fmt1,' saved successfully!'!itime/icheckpoint,'saved successfully!'
          ! write(*,*) 'Elapsed time (s)',real(trestart,4)
          ! write(*,*) 'Aproximated writing speed (MB/s)',real(((s3df*16.)*1e-6)/trestart,4)
          write(*,*) 'If necesseary restart from:',itime+1
       endif
    end if

  end subroutine restart
  
  subroutine init_restart_adios2()

    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use variables, only : numscalar
    use param, only : ilmn, nrhotime, ntime, mhd_active
    use var, only : itimescheme, iibm
    use mhd, only : mhd_equation
    
    implicit none

    integer :: ierror
    
    integer :: is
    character(len=80) :: varname
    
    call decomp_2d_init_io(io_restart)
    
    call decomp_2d_register_variable(io_restart, "ux", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart, "uy", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart, "uz", 1, 0, 0, mytype)

    call decomp_2d_register_variable(io_restart, "pp", 3, 0, 0, mytype, phG) !! XXX: need some way to handle the different grid here...

    do is = 1, numscalar
       write(varname,"('phi-',I2.2)") is
       call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
    end do

    if ((itimescheme.eq.2) .or. (itimescheme.eq.3)) then
       call decomp_2d_register_variable(io_restart, "dux-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duy-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duz-2", 1, 0, 0, mytype)

       do is = 1, numscalar
          write(varname,"('dphi-',I2.2,'-2')") is
          call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
       end do

       if (itimescheme.eq.3) then
          call decomp_2d_register_variable(io_restart, "dux-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duy-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duz-3", 1, 0, 0, mytype)

          do is = 1, numscalar
             write(varname,"('dphi-',I2.2,'-3')") is
             call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
          end do
       endif
    endif

    if (iibm .ne. 0) then
       call decomp_2d_register_variable(io_restart, "ep", 1, 0, 0, mytype)
    endif

    if (ilmn) then
       do is = 1, nrhotime
          write(varname,"('rho-',I2.2)") is
          call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
       end do
       do is = 1, ntime
          write(varname,"('drho-',I2.2)") is
          call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
       end do
       call decomp_2d_register_variable(io_restart, "mu", 1, 0, 0, mytype)
    end if
 
    if (mhd_active .and. mhd_equation == 'induction') then
       call decomp_2d_register_variable(io_restart, "bx", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "by", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "bz", 1, 0, 0, mytype)
       if (itimescheme.eq.3) then
          call decomp_2d_register_variable(io_restart, "dbx-2", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "dby-2", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "dbz-2", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "dbx-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "dby-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "dbz-3", 1, 0, 0, mytype)
       endif
    end if

  end subroutine init_restart_adios2
  !############################################################################
  !!  SUBROUTINE: apply_spatial_filter
  !############################################################################
  subroutine apply_spatial_filter(ux1,uy1,uz1,phi1)

    use param
    !use var, only: uxf1,uyf1,uzf1,uxf2,uyf2,uzf2,uxf3,uyf3,uzf3,di1,di2,di3,phif1,phif2,phif3
    use var, only: uxf1,uyf1,uzf1,di1,phif1
    use var, only: ux2,uy2,uz2, phi2,uxf2,uyf2,uzf2,di2,phif2
    use var, only: ux3,uy3,uz3, phi3,uxf3,uyf3,uzf3,di3,phif3
    use variables
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3), numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi11

    integer :: i,j,k,npaire

    !if (iscalar == 1) phi11=phi1(:,:,:,1) !currently only first scalar
    if (ifilter==1.or.ifilter==2) then
      call filx(uxf1,ux1,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0,ubcx)
      call filx(uyf1,uy1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
      call filx(uzf1,uz1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    else
      uxf1=ux1
      uyf1=uy1
      uzf1=uz1
      !if (iscalar == 1) phif1=phi11
    end if

    call transpose_x_to_y(uxf1,ux2)
    call transpose_x_to_y(uyf1,uy2)
    call transpose_x_to_y(uzf1,uz2)
    !if (iscalar == 1) call transpose_x_to_y(phif1,phi2)

    if (ifilter==1.or.ifilter==3) then ! all filter or y filter
      call fily(uxf2,ux2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx)
      call fily(uyf2,uy2,di2,fisy,fiffy,fifsy,fifwy,ysize(1),ysize(2),ysize(3),0,ubcy)
      call fily(uzf2,uz2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcz)
      !if (iscalar.eq.1) call fily(phif2,phi2,di2,fisy,fiffy,fifsy,fifwy,ysize(1),ysize(2),ysize(3),0)
    else
      uxf2=ux2
      uyf2=uy2
      uzf2=uz2
      !if (iscalar == 1) phif2=phi2
    end if

    call transpose_y_to_z(uxf2,ux3)
    call transpose_y_to_z(uyf2,uy3)
    call transpose_y_to_z(uzf2,uz3)
    !if (iscalar == 1) call transpose_y_to_z(phif2,phi3)

    if (ifilter==1.or.ifilter==2) then
      call filz(uxf3,ux3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
      call filz(uyf3,uy3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
      call filz(uzf3,uz3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0,ubcz)
      !if (iscalar.eq.1) call filz(phif3,phi3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0)
    else
      uxf3=ux3
      uyf3=uy3
      uzf3=uz3
      !if (iscalar == 1) phif3=phi3
    end if

    call transpose_z_to_y(uxf3,ux2)
    call transpose_z_to_y(uyf3,uy2)
    call transpose_z_to_y(uzf3,uz2)
    !if (iscalar == 1) call transpose_z_to_y(phif3,phi2)

    call transpose_y_to_x(ux2,ux1)
    call transpose_y_to_x(uy2,uy1)
    call transpose_y_to_x(uz2,uz1)
    !if (iscalar == 1) call transpose_y_to_x(phi2,phi11)

    !if (iscalar == 1) phi1(:,:,:,1)=phi11

  end subroutine apply_spatial_filter
  
  !############################################################################
  !!  SUBROUTINE: init_inflow_outflow
  !############################################################################
  subroutine init_inflow_outflow()

    use decomp_2d_io, only : decomp_2d_init_io, decomp_2d_register_variable

    use param, only : ntimesteps
     
    integer :: nplanes
    
    call decomp_2d_init_io(io_ioflow)

    nplanes = ntimesteps
    
    call decomp_2d_register_variable(io_ioflow, "ux", 1, 0, 1, mytype, opt_nplanes=nplanes)
    call decomp_2d_register_variable(io_ioflow, "uy", 1, 0, 1, mytype, opt_nplanes=nplanes)
    call decomp_2d_register_variable(io_ioflow, "uz", 1, 0, 1, mytype, opt_nplanes=nplanes)
    
  end subroutine init_inflow_outflow
  !############################################################################
  !!  SUBROUTINE: read_inflow
  !############################################################################
  subroutine read_inflow(ux1,uy1,uz1,ifileinflow)

    use decomp_2d_io
    use var, only: ux_inflow, uy_inflow, uz_inflow
    use param
    use MPI

    implicit none

    integer :: ifileinflow
    real(mytype), dimension(NTimeSteps,xsize(2),xsize(3)) :: ux1,uy1,uz1
    character(20) :: fninflow

    character(len=1024) :: inflow_file
    
    ! Recirculate inflows 
    if (ifileinflow>=ninflows) then 
      ifileinflow=mod(ifileinflow,ninflows)
    endif

    ! Read inflow
    write(fninflow,'(i20)') ifileinflow+1
    write(inflow_file, "(A)") trim(inflowpath)//'inflow'//trim(adjustl(fninflow))
    if (nrank==0) print *,'READING INFLOW FROM ',inflow_file
    
    call decomp_2d_open_io(io_ioflow, inflow_file, decomp_2d_read_mode)
    
    call decomp_2d_read_inflow(inflow_file,"ux",ntimesteps,ux_inflow,io_ioflow)
    call decomp_2d_read_inflow(inflow_file,"uy",ntimesteps,uy_inflow,io_ioflow)
    call decomp_2d_read_inflow(inflow_file,"uz",ntimesteps,uz_inflow,io_ioflow)

    call decomp_2d_close_io(io_ioflow, inflow_file)

  end subroutine read_inflow
  !############################################################################
  !!  SUBROUTINE: append_outflow
  !############################################################################
  subroutine append_outflow(ux,uy,uz,timestep)
 
    use decomp_2d_io
    use var, only: ux_recoutflow, uy_recoutflow, uz_recoutflow, ilist
    use param

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    integer, intent(in) :: timestep
    integer :: j,k

    if (nrank==0.and.mod(itime,ilist)==0) print *, 'Appending outflow', timestep 
    do k=1,xsize(3)
    do j=1,xsize(2)
      ux_recoutflow(timestep,j,k)=ux(xend(1),j,k)
      uy_recoutflow(timestep,j,k)=uy(xend(1),j,k)
      uz_recoutflow(timestep,j,k)=uz(xend(1),j,k)
    enddo
    enddo

    return
  end subroutine append_outflow
  !############################################################################
  !!  SUBROUTINE: write_outflow
  !############################################################################
  subroutine write_outflow(ifileoutflow)

    use decomp_2d_io
    use param
    use var, only: ux_recoutflow, uy_recoutflow, uz_recoutflow
    use MPI

    implicit none

    integer,intent(in) :: ifileoutflow
    character(20) :: fnoutflow
    character(80) :: outflow_file
    logical, save :: clean = .true.
    integer :: iomode

    logical :: dir_exists

    iomode = decomp_2d_write_mode
    
    write(fnoutflow,'(i20)') ifileoutflow
    write(outflow_file, "(A)") './out/inflow'//trim(adjustl(fnoutflow))
    if (nrank==0) print *,'WRITING OUTFLOW TO ', outflow_file
    
    call decomp_2d_open_io(io_ioflow, outflow_file, iomode)
    call decomp_2d_start_io(io_ioflow, outflow_file)

    call decomp_2d_write_outflow(outflow_file,"ux",ntimesteps,ux_recoutflow,io_ioflow)
    call decomp_2d_write_outflow(outflow_file,"uy",ntimesteps,uy_recoutflow,io_ioflow)
    call decomp_2d_write_outflow(outflow_file,"uz",ntimesteps,uz_recoutflow,io_ioflow)

    call decomp_2d_end_io(io_ioflow, outflow_file)
    call decomp_2d_close_io(io_ioflow, outflow_file)
    
  end subroutine write_outflow
  !############################################################################
  !##################################################################
  !!  SUBROUTINE: compute_cfldiff
  !! DESCRIPTION: Computes Diffusion/Fourier number
  !!      AUTHOR: Kay Schäfer
  !##################################################################
  subroutine compute_cfldiff()
     use param, only : xnu,dt,dx,dy,dz,istret
     use param, only : cfl_diff_sum, cfl_diff_x, cfl_diff_y, cfl_diff_z
     use param, only : mhd_active
     use variables, only : dyp
     use mhd, only: mhd_equation,rem

     implicit none

     cfl_diff_x = xnu * dt/ (dx**2)
     cfl_diff_z = xnu * dt/ (dz**2)

     if (istret == 0) then
        cfl_diff_y = xnu * dt / (dy**2)
     else
        cfl_diff_y = xnu * dt / (minval(dyp)**2)
     end if

     cfl_diff_sum = cfl_diff_x + cfl_diff_y + cfl_diff_z

     if (nrank==0) then
        write(*,*) '==========================================================='
        write(*,*) 'Diffusion number'
        write(*,"(' cfl_diff_x             :        ',F13.8)") cfl_diff_x
        write(*,"(' cfl_diff_y             :        ',F13.8)") cfl_diff_y
        write(*,"(' cfl_diff_z             :        ',F13.8)") cfl_diff_z
        write(*,"(' cfl_diff_sum           :        ',F13.8)") cfl_diff_sum
        write(*,*) '==========================================================='
     endif
     
     if( mhd_active .and. mhd_equation=='induction') then
 
        cfl_diff_x = dt/ (dx**2) / rem
        cfl_diff_z = dt/ (dz**2) / rem
   
        if (istret == 0) then
           cfl_diff_y = dt / (dy**2) / rem
        else
           cfl_diff_y = dt / (minval(dyp)**2) / rem
        end if
   
        cfl_diff_sum = cfl_diff_x + cfl_diff_y + cfl_diff_z
   
        if (nrank==0) then
           write(*,*) '==========================================================='
           write(*,*) 'Magnetic Diffusion number'
           write(*,"(' B cfl_diff_x             :        ',F13.8)") cfl_diff_x
           write(*,"(' B cfl_diff_y             :        ',F13.8)") cfl_diff_y
           write(*,"(' B cfl_diff_z             :        ',F13.8)") cfl_diff_z
           write(*,"(' B cfl_diff_sum           :        ',F13.8)") cfl_diff_sum
           write(*,*) '==========================================================='
        endif
     endif  

     return
  end subroutine compute_cfldiff
  !##################################################################
    !!  SUBROUTINE: compute_cfl
    !! DESCRIPTION: Computes CFl number for stretched mesh
    !!      AUTHOR: Kay Schäfer
  !##################################################################
  subroutine compute_cfl(ux,uy,uz)
    use param, only : dx,dy,dz,dt,istret
    use mpi
    use variables, only : dyp

    implicit none

    integer      :: code, i,j,k,jloc
    real(mytype) :: value_x, value_y, value_z, value_sum
    real(mytype) :: maxvalue_sum, maxvalue_sum_out, maxvalue_x, maxvalue_y,  maxvalue_z
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(4) :: cflmax_in, cflmax_out
    !
    maxvalue_x  =-1609._mytype
    maxvalue_y  =-1609._mytype
    maxvalue_z  =-1609._mytype
    maxvalue_sum=-1609._mytype
    !
    if (istret == 0) then
       do j = xstart(2), xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:)) / dx)
          value_y    = maxval(abs(uy(:,jloc,:)) / dy)
          value_z    = maxval(abs(uz(:,jloc,:)) / dz)
          value_sum  = maxval(abs(ux(:,jloc,:)) / dx + abs(uy(:,jloc,:)) / dy +    abs(uz(:,jloc,:)) / dz)
          !
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    else
       do j = xstart(2), xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:)) / dx)
          value_y    = maxval(abs(uy(:,jloc,:)) / dyp(j))
          value_z    = maxval(abs(uz(:,jloc,:)) / dz)
          value_sum  = maxval(abs(ux(:,jloc,:)) / dx + abs(uy(:,jloc,:)) / dyp(j) + abs(uz(:,jloc,:)) /dz)
          !
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    end if

    cflmax_in =  (/maxvalue_x, maxvalue_y, maxvalue_z, maxvalue_sum/)

    call    MPI_REDUCE(cflmax_in,cflmax_out,4,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    if (nrank == 0) then
      write(*,"(' CFL_x                  : ',F17.8)") cflmax_out(1) * dt
      write(*,"(' CFL_y                  : ',F17.8)") cflmax_out(2) * dt
      write(*,"(' CFL_z                  : ',F17.8)") cflmax_out(3) * dt
      !write(*,"(' CFL_sum                : ',F17.8)") cflmax_out(4)*dt
    end if
  end subroutine compute_cfl
  !##################################################################
  !##################################################################
  ! Rescale pressure to physical pressure
  ! Written by Kay Schäfer 2019
  !##################################################################
  elemental subroutine rescale_pressure(pre1)

    use param, only : itimescheme, gdt
    implicit none

    real(mytype), intent(inout) :: pre1

    ! Adjust pressure to physical pressure
    ! Multiply pressure by factor of time-scheme
    ! 1/gdt = 1  / (dt * c_k)
    !
    ! Explicit Euler, AB2, AB3, AB4, RK3
    if (itimescheme>=1 .and. itimescheme<=5) then
       pre1 = pre1 / gdt(3)
    ! RK4
    elseif (itimescheme==6) then
       pre1 = pre1 / gdt(5)
    endif

  end subroutine
  !##################################################################
  !##################################################################
  subroutine mean_plane_x (f1,nx,ny,nz,fm1)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f1
    real(mytype),intent(out),dimension(ny,nz) :: fm1
    integer :: i,j,k

    fm1 = sum(f1, DIM=1) / real(nx, mytype)
    return

  end subroutine mean_plane_x
  !##################################################################
  !##################################################################
  subroutine mean_plane_y (f2,nx,ny,nz,fm2)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f2
    real(mytype),intent(out),dimension(nx,nz) :: fm2
    integer :: i,j,k

    fm2 = sum(f2, DIM=2) / real(ny, mytype)
    return

  end subroutine mean_plane_y
  !##################################################################
  !##################################################################
  subroutine mean_plane_z (f3,nx,ny,nz,fm3)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f3
    real(mytype),intent(out),dimension(nx,ny) :: fm3
    integer :: i,j,k

    fm3 = sum(f3, DIM=3) / real(nz,mytype)
    return

  end subroutine mean_plane_z
  ! Subroutine to rename a file/directory
  !
  ! [string, in]            oldname  - name of existing file
  ! [string, in]            newname  - existing file to be renamed as
  ! [integer, in, optional] opt_rank - which rank should perform the operation? all others will wait.
  subroutine rename(oldname, newname, opt_rank)

    use MPI
    use decomp_2d_io, only : gen_iodir_name
    
    character(len=*), intent(in) :: oldname
    character(len=*), intent(in) :: newname
    integer, intent(in), optional :: opt_rank

    integer :: exe_rank
    logical :: exist
    character(len=:), allocatable :: cmd

    integer :: ierror

    character(len=:), allocatable :: oldname_ext
    
    if (present(opt_rank)) then
       exe_rank = opt_rank
    else
       exe_rank = 0
    end if

    if (nrank == exe_rank) then
       oldname_ext = gen_iodir_name(oldname, io_restart)
       inquire(file=oldname_ext, exist=exist)
       if (exist) then
          cmd = "mv "//oldname_ext//" "//newname
          call execute_command_line(cmd)
       end if
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    if (ierror /= 0) then
       call decomp_2d_abort(ierror, &
            "Error in MPI_Barrier!")
    end if
    
  end subroutine rename

  ! Subroutine to delete a file/director
  !
  ! [string, in]            name     - name of file to delete
  ! [integer, in, optional] opt_rank - which rank should perform the operation? all others will wait.
  subroutine delete_filedir(name, opt_rank)

    use MPI
    
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: opt_rank

    integer :: exe_rank
    logical :: exist
    integer :: unit

    character(len=:), allocatable :: cmd

    integer :: ierror
    
    if (present(opt_rank)) then
       exe_rank = opt_rank
    else
       exe_rank = 0
    end if

    if (nrank == exe_rank) then
       inquire(file=name, exist=exist)
       if (exist) then
          inquire(file=name//"/.", exist=exist)
          if (.not. exist) then
             ! Open file so it can be deleted on close
             open(newunit=unit, file=name)
             close(unit, status='delete')
          else
             ! Directory - hopefully rm should work...
             cmd = "rm -r "//name
             call execute_command_line(cmd)
          end if
       end if
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    if (ierror /= 0) then
       call decomp_2d_abort(ierror, &
            "Error in MPI_Barrier!")
    end if
    
  end subroutine delete_filedir

  ! Relatively crude validation function, checks that old and new checkpoints are the same size.
  ! XXX: for ADIOS2/BP4 this isn't really correct/useful as it's a directory...
  !
  ! [string, in]            refname  - name of previous checkpoint
  ! [string, in]            testname - name of new checkpoint
  ! [integer, in, optional] opt_rank - which rank should perform the operation? all others will wait.
  logical function validate_restart(refname, testname, opt_rank)

    use MPI
    use decomp_2d_io, only : gen_iodir_name
    
    character(len=*), intent(in) :: refname
    character(len=*), intent(in) :: testname
    integer, intent(in), optional :: opt_rank

    integer :: exe_rank
    logical :: success
    integer :: ierror

    integer :: refsize
    integer :: testsize

    logical :: refexist, testexist
    logical, save :: checked_initial = .false.

    character(len=:), allocatable :: testname_ext
    
    if (present(opt_rank)) then
       exe_rank = opt_rank
    else
       exe_rank = 0
    end if

    if (nrank == exe_rank) then
       success = .true.

       testname_ext = gen_iodir_name(testname, io_restart)
       
       inquire(file=refname, size=refsize, exist=refexist)
       inquire(file=testname_ext, size=testsize, exist=testexist)

       if (testexist) then
          if (refexist) then
             if (testsize /= refsize) then
                success = .false.
             end if
          else
             if (checked_initial) then
                print *, "ERROR: old restart ", refname, " doesn't exist!"
                success = .false.
             else
                ! Must assume this is the first call to restart, no old restart should exist
                success = .true.
             end if
             checked_initial = .true.
          end if
       else
          success = .false.
       end if
    end if

    call MPI_Bcast(success, 1, MPI_LOGICAL, exe_rank, MPI_COMM_WORLD, ierror)
    if (ierror /= 0) then
       call decomp_2d_abort(ierror, &
            "Error in MPI_Allreduce")
    end if

    validate_restart = success
    
  end function validate_restart
  
end module tools
!##################################################################

!===================================================
! Subroutine for computing the local and global CFL
! number, according to Lele 1992.
!===================================================
!##################################################################
subroutine cfl_compute(uxmax,uymax,uzmax)

  use decomp_2d_constants
  use param
  use variables
  use var

  implicit none

  real(mytype),intent(in) :: uxmax,uymax,uzmax
  real(mytype) :: cfl_x_adv,cfl_x_diff,cfl_y_adv,cfl_y_diff,cfl_z_adv,cfl_z_diff
  real(mytype) :: cfl_conv_lim, cfl_diff_lim
  real(mytype) :: sigma_conv(3), sigma_diff(3)
  real(mytype) :: visc

  ! Set the constants (this is true for periodic boundaries)
  sigma_conv=[zero, sqrt(three), 2.85_mytype]
  sigma_diff=[two, 2.5_mytype, 2.9_mytype]

  if(jles==0) then
     visc=xnu
  elseif (jles==1) then
     visc=xnu
  endif

  ! This is considering 1D peridic boundaries
  ! Do x-direction
  cfl_x_adv =abs(uxmax) * dt / dx
  cfl_x_diff = visc * dt / dx**2
  ! Do y-direction
  cfl_y_adv = abs(uymax) * dt / dy
  cfl_y_diff = visc * dt / dy**2
  ! Do z-direction
  cfl_z_adv = abs(uzmax) * dt / dz
  cfl_z_diff = visc * dt / dz**2

  ! So far we will focus on uniform grids
  if(nrank == 0) then
     write(*,*) ' '
     write(*,1002) cfl_x_adv, cfl_x_diff
1002 format('CFL x-direction (Adv and Diff) =',F9.4,',',F9.4)
     write(*,1003) cfl_y_adv, cfl_y_diff
1003 format('CFL y-direction (Adv and Diff) =',F9.4,',',F9.4)
     write(*,1004) cfl_z_adv, cfl_z_diff
1004 format('CFL z-direction (Adv and Diff) =',F9.4,',',F9.4)
     cfl_conv_lim = sigma_conv(itimescheme) / sqrt(three)
     cfl_diff_lim = sigma_diff(itimescheme) / six
     write(*,1005) cfl_conv_lim, cfl_diff_lim
     write(*,*) ' '
1005 format('CFL limits (Adv and Diff) : ',F9.4,',',F9.4)
  endif

end subroutine cfl_compute
!##################################################################
!##################################################################
subroutine inversion5_v1(aaa_in,eee,spI)

  use decomp_2d
  use decomp_2d_constants
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use mpi
  use utilities

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa, aaa_in
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  aaa = aaa_in

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, ny/2 - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k)=cx(tmp1,tmp2)
              eee(j,mi,k)=cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                             iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo

  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,ny/2,k,2)) / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,ny/2,k,2)) / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,ny/2,k,3)) - tmp1 * rl(aaa(j,ny/2-1,k,4)),&
                     iy(aaa(j,ny/2,k,3)) - tmp2 * iy(aaa(j,ny/2-1,k,4)))

        if (abs(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,ny/2,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,ny/2-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,ny/2,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,ny/2-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1,tmp2)
        eee(j,ny/2,k) = cx(tmp3,tmp4)

        if (abs(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1, tmp2)
        a1(j,k) = cx(rl(aaa(j,ny/2-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,ny/2-1,k,4)) * iy(b1(j,k)))
        eee(j,ny/2-1,k) = cx(rl(eee(j,ny/2-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,ny/2,k)),&
                             iy(eee(j,ny/2-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,ny/2,k)))
     enddo
  enddo

  do i = ny/2 - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one/iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) - rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) - iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v1
!##################################################################
!##################################################################
subroutine inversion5_v2(aaa,eee,spI)

  use decomp_2d
  use decomp_2d_constants
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use MPI
  use utilities

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, nym - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k) = cx(tmp1, tmp2)
              eee(j,mi,k) = cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                               iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo
  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,nym,k,2)) / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,nym,k,2)) / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,nym,k,3)) - tmp1 * rl(aaa(j,nym-1,k,4)),&
                     iy(aaa(j,nym,k,3)) - tmp2 * iy(aaa(j,nym-1,k,4)))
        if (abs(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,nym,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,nym-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,nym,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,nym-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1, tmp2)
        eee(j,nym,k) = cx(tmp3, tmp4)

        if (abs(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1,tmp2)
        a1(j,k) = cx(rl(aaa(j,nym-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,nym-1,k,4)) * iy(b1(j,k)))
        eee(j,nym-1,k) = cx(rl(eee(j,nym-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,nym,k)),&
                            iy(eee(j,nym-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,nym,k)))
     enddo
  enddo

  do i = nym - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one / iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) -rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) -iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v2
!##################################################################
!##################################################################
subroutine tripping(tb,ta)

  use param
  use variables
  use decomp_2d
  use decomp_2d_constants
  use decomp_2d_mpi
  use mpi

  implicit none

  integer :: i,j,k
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb, ta
  integer :: seed0, ii, code
  real(mytype) :: z_pos, randx, p_tr, b_tr, x_pos, y_pos!, A_tr

  !Done in X-Pencils
  seed0=randomseed !Seed for random number
  !A_tr=A_trip*min(1.0,0.8+real(itime)/200.0)
  !xs_tr=4.0/2.853
  !ys_tr=2.0/2.853
  !ts_tr=4.0/2.853
  !x0_tr=40.0/2.853
  A_tr = 0.1*dt

  if ((itime == ifirst).and.(nrank == 0)) then
     call random_seed(size=ii)
     call random_seed(put=seed0*(/ (1, i = 1, ii) /))

     !DEBUG:
     !call random_number(randx)
     !call MPI_BCAST(randx,1,real_type,0,MPI_COMM_WORLD,code)
     !write(*,*) 'RANDOM:', nrank, randx, ii
     !First random generation of h_nxt


     do j=1,z_modes

        call random_number(randx)
        h_coeff(j)=one*(randx-zpfive)
     enddo
     h_coeff=h_coeff/sqrt(real(z_modes,mytype))
  endif

  !Initialization h_nxt  (always bounded by xsize(3)^2 operations)
  if (itime == ifirst) then
     call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
     nxt_itr=0
     do k=1,xsize(3)
        h_nxt(k)=zero
        z_pos=-zlz*zpfive+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(two*pi*j*z_pos/zlz)
        enddo
     enddo
  end if



  !Time-loop
  i=int(t/ts_tr)
  if (i.ge.nxt_itr) then  !Nxt_itr is a global variable
     nxt_itr=i+1

     !First random generation of h
     h_i(:)=h_nxt(:)
     if (nrank  ==  0) then
        do j=1,z_modes
           call random_number(randx)
           h_coeff(j)=one*(randx-zpfive)
        enddo
        h_coeff=h_coeff/sqrt(real(z_modes,mytype)) !Non-dimensionalization
     end if

     call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)


     !Initialization h_nxt  (always bounded by z_steps^2 operations)
     do k=1,xsize(3)
        h_nxt(k)=zero
        z_pos=-zlz*zpfive+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(two*pi*j*z_pos/zlz)
        enddo
     enddo
  endif

  !Time coefficient
  p_tr=t/ts_tr-i
  b_tr=three*p_tr**2-two*p_tr**3

  !Creation of tripping velocity
  do i=1,xsize(1)
     x_pos=(xstart(1)+(i-1)-1)*dx
     do j=1,xsize(2)
        !y_pos=(xstart(2)+(j-1)-1)*dy
        y_pos=yp(xstart(2)+(j-1))
        do k=1,xsize(3)
           !g(z)*EXP_F(X,Y)
           ta(i,j,k)=((one-b_tr)*h_i(k)+b_tr*h_nxt(k))
           !ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-(y_pos/ys_tr)**2)*ta(i,j,k)
           ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-((y_pos-zpfive)/ys_tr)**2)*ta(i,j,k)
           tb(i,j,k)=tb(i,j,k)+ta(i,j,k)

           z_pos=-zlz*zpfive+(xstart(3)+(k-1)-1)*dz
           ! if ((((x_pos-x0_tr)**2).le.9.0e-3).and.(y_pos.le.0.0001).and.((z_pos).le.0.03))then
           !       open(442,file='tripping.dat',form='formatted',position='APPEND')
           !  write(442,*) t,ta(i,j,k)
           !  close(442)
           ! end if

        enddo
     enddo
  enddo

  return
end subroutine tripping
!##################################################################
!##################################################################
!!TRIPPING SUBROUTINE FOR TURBULENT BOUNDARY LAYERS
!##################################################################
subroutine tbl_tripping(tb,ta)

  use param
  use variables
  use decomp_2d
  use decomp_2d_constants
  use decomp_2d_mpi
  use mpi

  implicit none

  integer :: i,j,k
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb, ta
  integer :: seed0, ii, code
  !real(mytype) :: x0_tr_tbl, xs_tr_tbl,ys_tr_tbl,ts_tr_tbl !Scales related with maximum wave numbers
  real(mytype) :: z_pos, randx, p_tr,b_tr, x_pos, y_pos
  logical :: exist


  !Done in X-Pencils

  seed0=randomseed !Seed for random number
  !xs_tr_tbl=4.0/2.853
  !ys_tr_tbl=1.0/2.853
  !ts_tr_tbl=4.0/2.853
  !x0_tr_tbl=10.0/2.853

  !A_tr =  0.75/(ts_tr_tbl) !0.3/(ts_tr)


  if ((itime == ifirst).and.(nrank == 0)) then
     call random_seed(size=ii)
     call random_seed(put=seed0*(/ (1, i = 1, ii) /))

     INQUIRE(FILE='restart.nc',exist=exist)
     !if ((ilit==1).AND.(exist)) then
     !if (exist) then
     !write(*,*) 'h_coeff1 and phase1 already read from restart.nc'
     !write(*,*) 'h_coeff2 and phase2 already read from restart.nc'
     !nxt_itr=int(t/ts_tr_tbl)
     !else
     nxt_itr=1
     do j=1,z_modes
        call random_number(randx)
        h_coeff1(j)=one*(randx-zpfive)/sqrt(real(z_modes,mytype))
        call random_number(randx)
        phase1(j) = two*pi*randx
        call random_number(randx)
        h_coeff2(j)=one*(randx-zpfive)/sqrt(real(z_modes,mytype))
        call random_number(randx)
        phase2(j) = two*pi*randx
     enddo
     !endif
  endif

  !Initialization h_nxt  (always bounded by xsize(3)^2 operations)
  if (itime == ifirst) then
     call MPI_BCAST(h_coeff1,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(phase1,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(h_coeff2,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(phase2,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(nxt_itr,1,mpi_int,0,MPI_COMM_WORLD,code)

     do k=1,xsize(3)
        h_1(k)=zero
        h_2(k)=zero
        z_pos=-zlz*zpfive+real(xstart(3)+(k-1)-1,mytype)*dz
        do j=1,z_modes
           h_1(k)= h_1(k)+h_coeff1(j)*sin(two*pi*real(j,mytype)*z_pos/zlz+phase1(j))
           h_2(k)= h_2(k)+h_coeff2(j)*sin(two*pi*real(j,mytype)*z_pos/zlz+phase2(j))
        enddo
     enddo
  end if

  !Time-loop
  i=int(t/ts_tr_tbl)
  if (i.ge.nxt_itr) then
     nxt_itr=i+1
     !Move h_nxt to h_i
     h_2(:)=h_1(:)
     !---------------------------------------------------------
     !Create signal again
     if (nrank  ==  0) then
        do j=1,z_modes
           call random_number(randx)
           h_coeff1(j)=one*(randx-zpfive)/sqrt(real(z_modes,mytype))
           call random_number(randx)
           phase1(j) = two*pi*randx
        enddo
     end if

     call MPI_BCAST(h_coeff1,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(phase1,z_modes,real_type,0,MPI_COMM_WORLD,code)

     !Initialization h_nxt  (always bounded by z_steps^2 operations)
     do k=1,xsize(3)
        h_1(k)=zero
        z_pos=-zlz*zpfive+real(xstart(3)+(k-1)-1,mytype)*dz
        do j=1,z_modes
           h_1(k)= h_1(k)+h_coeff1(j)*sin(two*pi*real(j,mytype)*z_pos/zlz+phase1(j))
        enddo
     enddo
  endif
  !-------------------------------------------------------------------

  !Time coefficient
  p_tr=t/ts_tr_tbl-i
  b_tr=three*p_tr**2-two*p_tr**3
  !Creation of tripping velocity
  do i=1,xsize(1)
     x_pos=(xstart(1)+(i-1)-1)*dx
     do j=1,xsize(2)
        y_pos=yp(xstart(2)+(j-1))
        do k=1,xsize(3)
           ta(i,j,k)=((one-b_tr)*h_1(k)+b_tr*h_2(k))
           ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr_tbl)/xs_tr_tbl)**2-((y_pos-0.05_mytype)/ys_tr_tbl)**2)*ta(i,j,k)
           tb(i,j,k)=tb(i,j,k)+ta(i,j,k)

           z_pos=-zlz*zpfive+real(xstart(3)+(k-1)-1,mytype)*dz

        enddo
     enddo
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,code)
  !if (nrank==0) write(*,*) maxval(ta(:,:,:)),minval(ta), z_modes

  return
end subroutine tbl_tripping
!##################################################################
!##################################################################
subroutine calc_temp_eos(temp, rho, phi, mweight, xlen, ylen, zlen)

  use decomp_2d
  use decomp_2d_constants
  use param, only : pressure0, imultispecies
  use var, only : numscalar

  implicit none

  !! inputs
  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen) :: rho
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! outputs
  real(mytype), intent(out), dimension(xlen, ylen, zlen) :: temp

  !! locals
  real(mytype), dimension(xlen, ylen, zlen) :: mweight

  temp(:,:,:) = pressure0 / rho(:,:,:)
  if (imultispecies) then
     call calc_mweight(mweight, phi, xlen, ylen, zlen)
     temp(:,:,:) = temp(:,:,:) * mweight(:,:,:)
  endif

endsubroutine calc_temp_eos
!##################################################################
!##################################################################
subroutine calc_rho_eos(rho, temp, phi, mweight, xlen, ylen, zlen)

  use decomp_2d
  use decomp_2d_constants
  use param, only : pressure0, imultispecies
  use var, only : numscalar

  implicit none

  !! INPUTS
  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen) :: temp
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! OUTPUTS
  real(mytype), intent(out), dimension(xlen, ylen, zlen) :: rho

  !! LOCALS
  real(mytype), dimension(xlen, ylen, zlen) :: mweight

  rho(:,:,:) = pressure0 / temp(:,:,:)
  if (imultispecies) then
     call calc_mweight(mweight, phi, xlen, ylen, zlen)
     rho(:,:,:) = rho(:,:,:) * mweight(:,:,:)
  endif

endsubroutine calc_rho_eos
!##################################################################
!##################################################################
subroutine calc_mweight(mweight, phi, xlen, ylen, zlen)

  use decomp_2d
  use decomp_2d_constants
  use param, only : zero, one
  use param, only : massfrac, mol_weight
  use var, only : numscalar

  implicit none

  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! LOCALS
  real(mytype), dimension(xlen, ylen, zlen) :: mweight
  integer :: is

  mweight(:,:,:) = zero
  do is = 1, numscalar
     if (massfrac(is)) then
        mweight(:,:,:) = mweight(:,:,:) + phi(:,:,:,is) / mol_weight(is)
     endif
  enddo
  mweight(:,:,:) = one / mweight(:,:,:)

endsubroutine calc_mweight
!##################################################################
!##################################################################
subroutine test_min_max(name,text,array_tmp,i_size_array_tmp)

  use param
  use variables
  use decomp_2d
  use decomp_2d_constants
  use decomp_2d_mpi
  use MPI

  implicit none

  integer :: ierror, i, i_size_array_tmp
  real(mytype) :: max_tmp, min_tmp, tot_tmp
  real(mytype), dimension(i_size_array_tmp) :: array_tmp
  character(len=5) :: name
  character(len=15) :: text

  max_tmp=-0.000000000000000001_mytype
  tot_tmp=0._mytype
  min_tmp=+1000000000000000000._mytype
  do i=1,size(array_tmp)
    max_tmp=max(max_tmp,array_tmp(i))
    tot_tmp=tot_tmp + array_tmp(i)
    min_tmp=min(min_tmp,array_tmp(i))
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,max_tmp,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(MPI_IN_PLACE,min_tmp,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(MPI_IN_PLACE,tot_tmp,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  if (nrank == 0) then
     write(*,*) " "
     write(*,*) trim(text)//' Max ',name,max_tmp
     write(*,*) trim(text)//' Tot ',name,tot_tmp
     write(*,*) trim(text)//' Min ',name,min_tmp
     write(*,*) " "
     flush(6)
  endif

end subroutine test_min_max


