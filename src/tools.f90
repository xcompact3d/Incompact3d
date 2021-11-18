!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
module tools

  implicit none

  logical, save :: adios2_restart_initialised = .false.

  character(len=*), parameter :: io_restart = "restart-io"
  character(len=*), parameter :: resfile = "checkpoint"
  character(len=*), parameter :: io_ioflow = "in-outflow-io"
  
  private

  public :: test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, &
       apply_spatial_filter, read_inflow, append_outflow, write_outflow, init_inflow_outflow, &
       compute_cfldiff, compute_cfl, &
       rescale_pressure, mean_plane_x, mean_plane_y, mean_plane_z, &
       avg3d

contains
  !##################################################################
  !##################################################################
  subroutine test_scalar_min_max(phi)

    use decomp_2d
    use variables
    use param
    use var
    use mpi
    use dbg_schemes, only: abs_prec

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

        if (abs_prec(phimax1) > 100._mytype) then !if phi control turned off
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

    use decomp_2d
    use variables
    use param
    use var
    use mpi
    use dbg_schemes, only: abs_prec

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

       if((abs_prec(uxmax1)>=onehundred).or.(abs_prec(uymax1)>=onehundred).OR.(abs_prec(uzmax1)>=onehundred)) then
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

    use decomp_2d
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
          write(*,"(' Time step =',i7,'/',i7,', Time unit =',F9.4)") itime,ilast,t
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

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use navier, only : gradp

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
       call decomp_2d_write_one(3,pp3,resfile,"pp",0,io_restart,phG,reduce_prec=.false.)
       !
       if (iscalar==1) then
          do is=1, numscalar
             write(varname, *) "phi-", is
             call decomp_2d_write_one(1,phi1(:,:,:,is),resfile,varname,0,io_restart,reduce_prec=.false.)
             ! previous time-steps
             if ((itimescheme==2).or.(itimescheme==3)) then ! AB2 or AB3
                write(varname, *) "dphi-", is, "-2"
                call decomp_2d_write_one(1,dphi1(:,:,:,2,is),resfile,varname,0,io_restart,reduce_prec=.false.)
             end if
             !
             if (itimescheme==3) then ! AB3
               write(varname, *) "dphi-", is, "-3"
               call decomp_2d_write_one(1,dphi1(:,:,:,3,is),resfile,varname,0,io_restart,reduce_prec=.false.)
             end if
          end do
       endif
       if (ilmn) then
          do is = 1, nrhotime
             write(varname, *) "rho-", is
             call decomp_2d_write_one(1,rho1(:,:,:,is),resfile,varname,0,io_restart,reduce_prec=.false.)
          enddo
          do is = 1, ntime
             write(varname, *) "drho-", is
             call decomp_2d_write_one(1,drho1(:,:,:,is),resfile,varname,0,io_restart,reduce_prec=.false.)
          enddo
          call decomp_2d_write_one(1,mu1(:,:,:),resfile,"mu",0,io_restart,reduce_prec=.false.)
       endif

       call decomp_2d_end_io(io_restart, resfile)
       call decomp_2d_close_io(io_restart, resfile)

       ! Write info file for restart - Kay Schäfer
       if (nrank == 0) then
         write(filename,"('restart',I7.7,'.info')") itime
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
       call decomp_2d_read_one(3,pp3,resfile,"pp",io_restart,phG,reduce_prec=.false.)
       !
       if (iscalar==1) then
         do is=1, numscalar
            write(varname, *) "phi-", is
            call decomp_2d_read_one(1,phi1(:,:,:,is),resfile,varname,io_restart,reduce_prec=.false.)
           ! previous time-steps
           if ((itimescheme==2).or.(itimescheme==3)) then ! AB2 or AB3
             write(varname, *) "dphi-", is, "-2"
             call decomp_2d_read_one(1,dphi1(:,:,:,2,is),resfile,varname,io_restart,reduce_prec=.false.)
           end if
           !
           if (itimescheme==3) then ! AB3
              write(varname, *) "dphi-", is, "-3"
              call decomp_2d_read_one(1,dphi1(:,:,:,3,is),resfile,varname,io_restart,reduce_prec=.false.)
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
             write(varname, *) "rho-", is
             call decomp_2d_read_one(1,rho1(:,:,:,is),resfile,varname,io_restart,reduce_prec=.false.)
          enddo
          do is = 1, ntime
             write(varname, *) "drho-", is
             call decomp_2d_read_one(1,drho1(:,:,:,is),resfile,varname,io_restart,reduce_prec=.false.)
          enddo
          call decomp_2d_read_one(1,mu1,resfile,"mu",io_restart,reduce_prec=.false.)
       end if

       call decomp_2d_end_io(io_restart, resfile)
       call decomp_2d_close_io(io_restart, resfile)

       !! Read time of restart file
       write(filename,"('restart',I7.7,'.info')") ifirst-1
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
         itime0 = ifirst-1
       end if
       
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

    use decomp_2d, only : mytype, phG
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use variables, only : numscalar
    use param, only : ilmn, nrhotime, ntime
    use var, only : itimescheme, iibm
    
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
       write(varname,*) "phi-", is
       call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
    end do

    if ((itimescheme.eq.2) .or. (itimescheme.eq.3)) then
       call decomp_2d_register_variable(io_restart, "dux-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duy-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duz-2", 1, 0, 0, mytype)

       do is = 1, numscalar
          write(varname,*) "dphi-", is, "-2"
          call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
       end do

       if (itimescheme.eq.3) then
          call decomp_2d_register_variable(io_restart, "dux-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duy-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duz-3", 1, 0, 0, mytype)

          do is = 1, numscalar
             write(varname,*) "dphi-", is, "-3"
             call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
          end do
       endif
    endif

    if (iibm .ne. 0) then
       call decomp_2d_register_variable(io_restart, "ep", 1, 0, 0, mytype)
    endif

    if (ilmn) then
       do is = 1, nrhotime
          write(varname, *) "rho-", is
          call decomp_2d_register_variable(io_restart, varname, 1, 0, 0, mytype)
       end do
       do is = 1, ntime
          write(varname, *) "drho-", is
          call decomp_2d_register_variable(io_restart, varname, 1, 0, 0, mytype)
       end do
    end if
    
  end subroutine init_restart_adios2
  !############################################################################
  !!  SUBROUTINE: apply_spatial_filter
  !############################################################################
  subroutine apply_spatial_filter(ux1,uy1,uz1,phi1)

    use decomp_2d
    use param
    use var, only: uxf1,uyf1,uzf1,uxf2,uyf2,uzf2,uxf3,uyf3,uzf3,di1,di2,di3,phif1,phif2,phif3
    use variables
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3), numscalar), intent(inout) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi11
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2, phi2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3, phi3

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

    use decomp_2d, only : mytype
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

    use decomp_2d
    use decomp_2d_io
    use var, only: ux_inflow, uy_inflow, uz_inflow
    use param
    use MPI

    implicit none

    integer :: ifileinflow
    real(mytype), dimension(NTimeSteps,xsize(2),xsize(3)) :: ux1,uy1,uz1
    character(20) :: fninflow

    character(80) :: inflow_file
    
    ! Recirculate inflows 
    if (ifileinflow>=ninflows) then 
      ifileinflow=mod(ifileinflow,ninflows)
    endif

    ! Read inflow
    write(fninflow,'(i20)') ifileinflow+1
#ifndef ADIOS2
    write(inflow_file, "(A)") trim(inflowpath)//'inflow'//trim(adjustl(fninflow))
#else
    write(inflow_file, "(A)") trim(inflowpath)//'inflow'
#endif
    if (nrank==0) print *,'READING INFLOW FROM ',inflow_file
    
    call decomp_2d_open_io(io_ioflow, inflow_file, decomp_2d_read_mode)
    call decomp_2d_start_io(io_ioflow, inflow_file)

    call decomp_2d_read_inflow(inflow_file,"ux",ntimesteps,ux_inflow,io_ioflow)
    call decomp_2d_read_inflow(inflow_file,"uy",ntimesteps,uy_inflow,io_ioflow)
    call decomp_2d_read_inflow(inflow_file,"uz",ntimesteps,uz_inflow,io_ioflow)

    call decomp_2d_end_io(io_ioflow, inflow_file)
    call decomp_2d_close_io(io_ioflow, inflow_file)

  end subroutine read_inflow
  !############################################################################
  !!  SUBROUTINE: append_outflow
  !############################################################################
  subroutine append_outflow(ux,uy,uz,timestep)
 
    use decomp_2d
    use decomp_2d_io
    use var, only: ux_recoutflow, uy_recoutflow, uz_recoutflow
    use param

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    integer, intent(in) :: timestep
    integer :: j,k

    if (nrank==0) print *, 'Appending outflow', timestep 
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

    use decomp_2d
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

    if (clean .and. (irestart .eq. 0)) then
       iomode = decomp_2d_write_mode
       clean = .false.
    else
       iomode = decomp_2d_append_mode
    end if
    
    write(fnoutflow,'(i20)') ifileoutflow
#ifndef ADIOS2
    write(outflow_file, "(A)") './out/inflow'//trim(adjustl(fnoutflow))
#else
    write(outflow_file, "(A)") './out/inflow'
#endif
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
     use variables, only : dyp
     use decomp_2d, only : nrank

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

     return
  end subroutine compute_cfldiff
  !##################################################################
    !!  SUBROUTINE: compute_cfl
    !! DESCRIPTION: Computes CFl number for stretched mesh
    !!      AUTHOR: Kay Schäfer
  !##################################################################
  subroutine compute_cfl(ux,uy,uz)
    use param, only : dx,dy,dz,dt,istret
    use decomp_2d, only : nrank, mytype, xsize, xstart, xend, real_type
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

    use decomp_2d, only : mytype
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
  !############################################################################
  !!
  !!  SUBROUTINE: avg3d
  !!      AUTHOR: Stefano Rolfo
  !! DESCRIPTION: Compute the total sum of a a 3d field
  !!
  !############################################################################
  subroutine avg3d (var, avg)

    use decomp_2d, only: real_type, xsize, xend
    use param
    use variables, only: nx,ny,nz,nxm,nym,nzm
    use mpi

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: var
    real(mytype), intent(out) :: avg
    real(mytype)              :: dep

    integer :: i,j,k, code
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

    if (nclx1==1.and.xend(1)==nx) then
       xsize1=xsize(1)-1
    else
       xsize1=xsize(1)
    endif
    if (ncly1==1.and.xend(2)==ny) then
       xsize2=xsize(2)-1
    else
       xsize2=xsize(2)
    endif
    if (nclz1==1.and.xend(3)==nz) then
       xsize3=xsize(3)-1
    else
       xsize3=xsize(3)
    endif
    if (nclx1==1) then
       nxc=nxm
    else
       nxc=nx
    endif
    if (ncly1==1) then
       nyc=nym
    else
       nyc=ny
    endif
    if (nclz1==1) then
       nzc=nzm
    else
       nzc=nz
    endif

    dep=zero
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             !dep=dep+var(i,j,k)**2
             dep=dep+var(i,j,k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(dep,avg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    avg=avg/(nxc*nyc*nzc)

    return

  end subroutine avg3d
end module tools
!##################################################################

!===================================================
! Subroutine for computing the local and global CFL
! number, according to Lele 1992.
!===================================================
!##################################################################
subroutine cfl_compute(uxmax,uymax,uzmax)

  use param
  use variables
  use var
  use dbg_schemes, only: sqrt_prec, abs_prec

  implicit none

  real(mytype),intent(in) :: uxmax,uymax,uzmax
  real(mytype) :: cfl_x_adv,cfl_x_diff,cfl_y_adv,cfl_y_diff,cfl_z_adv,cfl_z_diff
  real(mytype) :: cfl_conv_lim, cfl_diff_lim
  real(mytype) :: sigma_conv(3), sigma_diff(3)
  real(mytype) :: visc

  ! Set the constants (this is true for periodic boundaries)
  sigma_conv=[zero, sqrt_prec(three), 2.85_mytype]
  sigma_diff=[two, 2.5_mytype, 2.9_mytype]

  if(jles==0) then
     visc=xnu
  elseif (jles==1) then
     visc=xnu
  endif

  ! This is considering 1D peridic boundaries
  ! Do x-direction
  cfl_x_adv =abs_prec(uxmax) * dt / dx
  cfl_x_diff = visc * dt / dx**2
  ! Do y-direction
  cfl_y_adv = abs_prec(uymax) * dt / dy
  cfl_y_diff = visc * dt / dy**2
  ! Do z-direction
  cfl_z_adv = abs_prec(uzmax) * dt / dz
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
     cfl_conv_lim = sigma_conv(itimescheme) / sqrt_prec(three)
     cfl_diff_lim = sigma_diff(itimescheme) / six
     write(*,1005) cfl_conv_lim, cfl_diff_lim
     write(*,*) ' '
1005 format('CFL limits (Adv and Diff) : ',F9.4,',',F9.4)
  endif

end subroutine cfl_compute
!##################################################################
!##################################################################
subroutine stretching()

  use decomp_2d
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use mpi
  use dbg_schemes, only: abs_prec, sqrt_prec, sin_prec, cos_prec, tan_prec, atan_prec

  implicit none

  real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
  integer :: j

  yinf=-yly/two
  den=two*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)
  xcx=one/beta/alpha
  if (alpha.ne.0.) then
     if (istret.eq.1) yp(1)=zero
     if (istret.eq.2) yp(1)=zero
     if (istret.eq.1) yeta(1)=zero
     if (istret.eq.2) yeta(1)=-half
     if (istret.eq.3) yp(1)=zero
     if (istret.eq.3) yeta(1)=-half
     do j=2,ny
        if (istret==1) yeta(j)=real(j-1,mytype)*(one/nym)
        if (istret==2) yeta(j)=real(j-1,mytype)*(one/nym)-half
        if (istret==3) yeta(j)=real(j-1,mytype)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=two*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yeta(j))+one
        xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst-yinf
           if (yeta(j).eq.half) yp(j)=zero-yinf
           if (yeta(j).gt.half) yp(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst+yly
           if (yeta(j).eq.half) yp(j)=zero+yly
           if (yeta(j).gt.half) yp(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yeta(j).lt.half) yp(j)=(xnum1-cst+yly)*two
           if (yeta(j).eq.half) yp(j)=(zero+yly)*two
           if (yeta(j).gt.half) yp(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     yp(1)=-1.e10
     do j=2,ny
        yeta(j)=real(j-1,mytype)*(one/ny)
        yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
     enddo
  endif
  if (alpha.ne.0.) then
     do j=1,ny
        if (istret==1) yetai(j)=(real(j,mytype)-half)*(one/nym)
        if (istret==2) yetai(j)=(real(j,mytype)-half)*(one/nym)-half
        if (istret==3) yetai(j)=(real(j,mytype)-half)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yetai(j))+one
        xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst-yinf
           if (yetai(j).eq.half) ypi(j)=zero-yinf
           if (yetai(j).gt.half) ypi(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst+yly
           if (yetai(j).eq.half) ypi(j)=zero+yly
           if (yetai(j).gt.half) ypi(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yetai(j).lt.half) ypi(j)=(xnum1-cst+yly)*two
           if (yetai(j).eq.half) ypi(j)=(zero+yly)*two
           if (yetai(j).gt.half) ypi(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     ypi(1)=-1.e10
     do j=2,ny
        yetai(j)=real(j-1,mytype)*(one/ny)
        ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
     enddo
  endif

  !Mapping!!, metric terms
  if (istret .ne. 3) then
     do j=1,ny
        ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
        pp2y(j)=ppy(j)*ppy(j)
        pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
     enddo
     do j=1,ny
        ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
        pp2yi(j)=ppyi(j)*ppyi(j)
        pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
     enddo
  endif

  if (istret .eq. 3) then
     do j=1,ny
        ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
        pp2y(j)=ppy(j)*ppy(j)
        pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))/two
     enddo
     do j=1,ny
        ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
        pp2yi(j)=ppyi(j)*ppyi(j)
        pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))/two
     enddo
  endif

  !   yp(1) = 0.0
  !   yp(2) = 0.01
  !   coeff0= 1.1
  !   blender1 = 0.0
  !   blender2 = 0.0
  !   do j=3,ny
  !!      yeta(j)=(j-1.)*(1./ny)
  !!      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
  !
  !     if (yp(j-1).LE.3.5*1.0) then
  !       dy_plus_target = 8.0
  !       !Calculate re_tau guess somewhere
  !      dy_plus_current= (yp(j-1)-yp(j-2))*85.0
  !       !dy_plus_coeff is from 1 to 0
  !       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
  !       coeff = coeff0**dy_plus_coeff
  !
  !       dy_plus_coeff_old1 = dy_plus_coeff   !will be required for blenders
  !     else if (yp(j-1).GE.39.0*1.0) then
  !       dy_plus_target = 10.0
  !       !Calculate re_tau guess somewhere
  !       dy_plus_current= (yp(j-1)-yp(j-2))*85.0
  !       !dy_plus_coeff is from 1 to 0
  !       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
  !
  !       if (blender2.LT.1.0) blender2 = blender2 + 0.1   !carry the coeff smoothly
  !       coeff = coeff0**((1.0-blender2)*dy_plus_coeff_old2+blender2*dy_plus_coeff)
  !     else
  !       dy_plus_target = 80.0
  !       !Calculate re_tau guess somewhere
  !       dy_plus_current= (yp(j-1)-yp(j-2))*85.0
  !       !dy_plus_coeff is from 1 to 0
  !       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
  !
  !       if (blender1.LT.1.0) blender1 = blender1 + 0.1   !carry the coeff smoothly
  !       coeff = coeff0**((1.0-blender1)*dy_plus_coeff_old1+blender1*dy_plus_coeff)
  !
  !       dy_plus_coeff_old2 = dy_plus_coeff   !will be required for blenders
  !     endif
  !     yp(j) = yp(j-1)+(yp(j-1)-yp(j-2))*coeff
  !   enddo
  !
  !   !Normalize to yly
  !   ypmax = yp(ny)
  !   yp = yp/ypmax*yly

  if (nrank == 0) then
     open(10,file='yp.dat', form='formatted')
     do j=1,ny
        write(10,*)yp(j)
     enddo
     close(10)
     open(10,file='ypi.dat', form='formatted')
     do j=1,nym
        write(10,*)ypi(j)
     enddo
     close(10)
  endif

end subroutine stretching
!##################################################################
!##################################################################
subroutine inversion5_v1(aaa_in,eee,spI)

  use decomp_2d
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use mpi
  use dbg_schemes, only: abs_prec

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

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

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
        if (abs_prec(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,ny/2,k,2)) / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,ny/2,k,2)) / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,ny/2,k,3)) - tmp1 * rl(aaa(j,ny/2-1,k,4)),&
                     iy(aaa(j,ny/2,k,3)) - tmp2 * iy(aaa(j,ny/2-1,k,4)))

        if (abs_prec(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,ny/2,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,ny/2-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs_prec(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,ny/2,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,ny/2-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1,tmp2)
        eee(j,ny/2,k) = cx(tmp3,tmp4)

        if (abs_prec(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
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
           if (abs_prec(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs_prec(iy(aaa(j,i,k,3))) > epsilon) then
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
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use MPI
  use dbg_schemes, only: abs_prec

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

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

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
        if (abs_prec(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,nym,k,2)) / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,nym,k,2)) / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,nym,k,3)) - tmp1 * rl(aaa(j,nym-1,k,4)),&
                     iy(aaa(j,nym,k,3)) - tmp2 * iy(aaa(j,nym-1,k,4)))
        if (abs_prec(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,nym,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,nym-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs_prec(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,nym,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,nym-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1, tmp2)
        eee(j,nym,k) = cx(tmp3, tmp4)

        if (abs_prec(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,nym-1,k,3))) > epsilon) then
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
           if (abs_prec(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs_prec(iy(aaa(j,i,k,3))) > epsilon) then
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
  use mpi
  use dbg_schemes, only: sqrt_prec, sin_prec, exp_prec

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
     h_coeff=h_coeff/sqrt_prec(real(z_modes,mytype))
  endif

  !Initialization h_nxt  (always bounded by xsize(3)^2 operations)
  if (itime == ifirst) then
     call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
     nxt_itr=0
     do k=1,xsize(3)
        h_nxt(k)=zero
        z_pos=-zlz*zpfive+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_nxt(k)= h_nxt(k)+h_coeff(j)*sin_prec(two*pi*j*z_pos/zlz)
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
        h_coeff=h_coeff/sqrt_prec(real(z_modes,mytype)) !Non-dimensionalization
     end if

     call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)


     !Initialization h_nxt  (always bounded by z_steps^2 operations)
     do k=1,xsize(3)
        h_nxt(k)=zero
        z_pos=-zlz*zpfive+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_nxt(k)= h_nxt(k)+h_coeff(j)*sin_prec(two*pi*j*z_pos/zlz)
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
           !ta(i,j,k)=A_tr*exp_prec(-((x_pos-x0_tr)/xs_tr)**2-(y_pos/ys_tr)**2)*ta(i,j,k)
           ta(i,j,k)=A_tr*exp_prec(-((x_pos-x0_tr)/xs_tr)**2-((y_pos-zpfive)/ys_tr)**2)*ta(i,j,k)
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
  use mpi
  use dbg_schemes, only: sqrt_prec, exp_prec, sin_prec

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
        h_coeff1(j)=one*(randx-zpfive)/sqrt_prec(real(z_modes,mytype))
        call random_number(randx)
        phase1(j) = two*pi*randx
        call random_number(randx)
        h_coeff2(j)=one*(randx-zpfive)/sqrt_prec(real(z_modes,mytype))
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
           h_1(k)= h_1(k)+h_coeff1(j)*sin_prec(two*pi*real(j,mytype)*z_pos/zlz+phase1(j))
           h_2(k)= h_2(k)+h_coeff2(j)*sin_prec(two*pi*real(j,mytype)*z_pos/zlz+phase2(j))
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
           h_coeff1(j)=one*(randx-zpfive)/sqrt_prec(real(z_modes,mytype))
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
           h_1(k)= h_1(k)+h_coeff1(j)*sin_prec(two*pi*real(j,mytype)*z_pos/zlz+phase1(j))
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
           ta(i,j,k)=A_tr*exp_prec(-((x_pos-x0_tr_tbl)/xs_tr_tbl)**2-((y_pos-0.05_mytype)/ys_tr_tbl)**2)*ta(i,j,k)
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
function rl(complexnumber)

  use param

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!##################################################################
!##################################################################
function iy(complexnumber)

  use param

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!##################################################################
!##################################################################
function cx(realpart,imaginarypart)

  use param

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!##################################################################
!##################################################################
subroutine calc_temp_eos(temp, rho, phi, mweight, xlen, ylen, zlen)

  use decomp_2d
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
function r8_random ( s1, s2, s3 )

!*****************************************************************************80
!
!! R8_RANDOM returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function returns a pseudo-random number rectangularly distributed
!    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
!    of Applied Statistics (1984) volume 33), not as claimed in the
!    original article.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    FORTRAN77 original version by Brian Wichman, David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Brian Wichman, David Hill,
!    Algorithm AS 183: An Efficient and Portable Pseudo-Random
!    Number Generator,
!    Applied Statistics,
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
!    seed for the sequence.  These values should be positive
!    integers between 1 and 30,000.
!
!    Output, real ( kind = 8 ) R8_RANDOM, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) s3
  real ( kind = 8 ) r8_random

  s1 = mod ( 171 * s1, 30269 )
  s2 = mod ( 172 * s2, 30307 )
  s3 = mod ( 170 * s3, 30323 )

  r8_random = mod ( real ( s1, kind = 8 ) / 30269.0D+00 &
                  + real ( s2, kind = 8 ) / 30307.0D+00 &
                  + real ( s3, kind = 8 ) / 30323.0D+00, 1.0D+00 )

  return
end
!##################################################################
function return_30k(x) result(y)

  integer ( kind = 4 ), intent(in) :: x
  integer ( kind = 4 )             :: y
  integer ( kind = 4 ), parameter  :: xmax = 30000

  y = iabs(x) - int(iabs(x)/xmax)*xmax
end function return_30k
!##################################################################
function r8_uni ( s1, s2 )

!*****************************************************************************80
!
!! R8_UNI returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function generates uniformly distributed pseudorandom numbers
!    between 0 and 1, using the 32-bit generator from figure 3 of
!    the article by L'Ecuyer.
!
!    The cycle length is claimed to be 2.30584E+18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original Pascal original version by Pierre L'Ecuyer
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Pierre LEcuyer,
!    Efficient and Portable Combined Random Number Generators,
!    Communications of the ACM,
!    Volume 31, Number 6, June 1988, pages 742-751.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
!    seed for the sequence.  On first call, the user should initialize
!    S1 to a value between 1 and 2147483562;  S2 should be initialized
!    to a value between 1 and 2147483398.
!
!    Output, real ( kind = 8 ) R8_UNI, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uni
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) z

  k = s1 / 53668
  s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
  if ( s1 < 0 ) then
    s1 = s1 + 2147483563
  end if

  k = s2 / 52774
  s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
  if ( s2 < 0 ) then
    s2 = s2 + 2147483399
  end if

  z = s1 - s2
  if ( z < 1 ) then
    z = z + 2147483562
  end if

  r8_uni = real ( z, kind = 8 ) / 2147483563.0D+00

  return
end
!##################################################################
!##################################################################
subroutine test_min_max(name,text,array_tmp,i_size_array_tmp)

  use param
  use variables
  use decomp_2d
  use MPI

  implicit none

  integer :: ierror, i, i_size_array_tmp
  real(mytype) :: max_tmp, min_tmp, tot_tmp, max_tot, min_tot, tot_tot
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
  call MPI_ALLREDUCE(max_tmp,max_tot,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(min_tmp,min_tot,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(tot_tmp,tot_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  if (nrank == 0) then
     write(*,*) " "
     write(*,*) trim(text)//' Max ',name,max_tot
     write(*,*) trim(text)//' Tot ',name,tot_tot
     write(*,*) trim(text)//' Min ',name,min_tot
     write(*,*) " "
     call flush(6)
  endif

  return
end subroutine test_min_max
