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

  private

  public :: test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, &
       compute_cfldiff, compute_cfl, &
       rescale_pressure, mean_plane_x, mean_plane_y, mean_plane_z, &
       channel_cfr

contains
  !##################################################################
  !##################################################################
  subroutine test_scalar_min_max(phi)

    USE decomp_2d
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    integer :: code,ierror,i,j,k,is,jglob
    real(mytype) :: phimax,phimin,phimax1,phimin1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(2,numscalar) :: phimaxin,phimaxout

    do is=1, numscalar

      ta1(:,:,:) = phi(:,:,:,is)
      ! ibm
      if (iibm.gt.0) then
        ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
      endif

      phimax=-1609.; phimin=1609.
      phimax = maxval(ta1(:,:,:))
      phimin =-minval(ta1(:,:,:))
      phimaxin(:,is) =  (/phimin, phimax /)
    enddo

    !call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    !call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(phimaxin,phimaxout,numscalar*2,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    do is=1,numscalar
      if (nrank.eq.0) then
        phimin1 = -phimaxout(1,is)
        phimax1 =  phimaxout(2,is)

        print *,'Phi'//char(48+is)//' min max=', real(phimin1,4), real(phimax1,4)

        if (abs(phimax1).ge.10.) then !if phi control turned off
           print *,'Scalar diverged! SIMULATION IS STOPPED!'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
        endif
      endif

    enddo

    return
  end subroutine test_scalar_min_max
  !##################################################################
  !##################################################################
  subroutine test_speed_min_max(ux,uy,uz)

    USE decomp_2d
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    integer :: code,ierror,i,j,k
    real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin
    real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(6) :: umaxin, umaxout

    if (iibm.gt.0) then
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

    if (nrank.eq.0) then

       print *,'U,V,W min=',real(uxmin1,4),real(uymin1,4),real(uzmin1,4)
       print *,'U,V,W max=',real(uxmax1,4),real(uymax1,4),real(uzmax1,4)
       !print *,'CFL=',real(abs(max(uxmax1,uymax1,uzmax1)*dt)/min(dx,dy,dz),4)

       if((abs(uxmax1).ge.10.).OR.(abs(uymax1).ge.10.).OR.(abs(uzmax1).ge.10.)) then
         print *,'Velocity diverged! SIMULATION IS STOPPED!'
         call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif

    endif

    return
  end subroutine test_speed_min_max
  !##################################################################
  !##################################################################
  subroutine simu_stats(iwhen)

    USE decomp_2d
    USE simulation_stats
    USE var
    USE MPI

    implicit none

    integer :: iwhen

    if (iwhen.eq.1) then !AT THE START OF THE SIMULATION
       tstart=zero;time1=zero;trank=zero;tranksum=zero;ttotal=zero
       call cpu_time(tstart)
    else if (iwhen.eq.2) then !AT THE START OF A TIME STEP
       call cpu_time(time1)
       if (nrank==0) then
          print *,'==========================================================='
          write(*,"(' Time step =',i7,'/',i7,', Time unit =',F9.4)") itime,ilast,t
       endif
    else if ((iwhen.eq.3).and.(itime.gt.ifirst)) then !AT THE END OF A TIME STEP
       call cpu_time(trank)
       if (nrank==0) print *,'Time for this time step (s):',real(trank-time1)
       telapsed = (trank-tstart)/thirtysixthousand
       tremaining  = telapsed*(ilast-itime)/(itime-ifirst)
       if (nrank==0) then
          write(*,"(' Remaining time:',I8,' h ',I2,' min')") int(tremaining), int((tremaining-int(tremaining))*sixty)
          write(*,"(' Elapsed time:  ',I8,' h ',I2,' min')") int(telapsed), int((telapsed-int(telapsed))*sixty)
       endif
    else if (iwhen.eq.4) then !AT THE END OF THE SIMULATION
       call cpu_time(trank); ttotal=trank-tstart
       if (nrank==0) then
          print *,'==========================================================='
          print *,''
          print *,'Good job! Xcompact3d finished successfully!'
          print *,''
          print *,'2DECOMP with p_row*p_col=',p_row,p_col
          print *,''
          print *,'nx*ny*nz=',nx*ny*nz
          print *,'nx,ny,nz=',nx,ny,nz
          print *,'dx,dy,dz=',dx,dy,dz
          print *,''
          print *,'Averaged time per step (s):',real(ttotal/(ilast-(ifirst-1)),4)
          print *,'Total wallclock (s):',real(ttotal,4)
          print *,'Total wallclock (m):',real(ttotal/sixty,4)
          print *,'Total wallclock (h):',real(ttotal/thirtysixthousand,4)
          print *,''
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
  subroutine restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,dphi1,px1,py1,pz1,iresflg)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
    use navier, only : gradp

    implicit none

    integer :: i,j,k,iresflg,nzmsize,fh,ierror,is,it,code
    integer :: ierror_o=0 !error to open sauve file during restart
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    real(mytype) :: xdt,tfield
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods
    logical :: fexists
    character(len=30) :: filename, filestart
    character(len=32) :: fmt2,fmt3,fmt4
    character(len=7) :: fmt1
    NAMELIST /Time/ tfield, itime
    NAMELIST /NumParam/ nx, ny, nz, istret, beta, dt, itimescheme

    write(filename,"('restart',I7.7)") itime
    write(filestart,"('restart',I7.7)") ifirst-1

    if (iresflg .eq. 1 ) then !Writing restart
       if (mod(itime, icheckpoint).ne.0) then
          return
       endif

       if (nrank==0) then
          print *,'===========================================================<<<<<'
          print *,'Writing restart point ',filename !itime/icheckpoint
          ! print *,'File size',real((s3df*16.)*1e-9,4),'GB'
       endif
    end if

    if (iresflg==1) then !write
       call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call decomp_2d_write_var(fh,disp,1,ux1)
       call decomp_2d_write_var(fh,disp,1,uy1)
       call decomp_2d_write_var(fh,disp,1,uz1)
       ! write previous time-step if necessary for AB2 or AB3
       if ((itimescheme.eq.2).or.(itimescheme.eq.3).or.(itimescheme.eq.7)) then
         call decomp_2d_write_var(fh,disp,1,dux1(:,:,:,2))
         call decomp_2d_write_var(fh,disp,1,duy1(:,:,:,2))
         call decomp_2d_write_var(fh,disp,1,duz1(:,:,:,2))
       end if
       ! for AB3 one more previous time-step
       if ((itimescheme.eq.3).or.(itimescheme.eq.7)) then
         call decomp_2d_write_var(fh,disp,1,dux1(:,:,:,3))
         call decomp_2d_write_var(fh,disp,1,duy1(:,:,:,3))
         call decomp_2d_write_var(fh,disp,1,duz1(:,:,:,3))
       end if
       !
       call decomp_2d_write_var(fh,disp,3,pp3,phG)
       !
       if (iscalar==1) then
          do is=1, numscalar
             call decomp_2d_write_var(fh,disp,1,phi1(:,:,:,is))
             ! previous time-steps
             if ((itimescheme.eq.2).or.(itimescheme.eq.3).or.(itimescheme.eq.7)) then ! AB2 or AB3
               call decomp_2d_write_var(fh,disp,1,dphi1(:,:,:,2,is))
             end if
             !
             if ((itimescheme.eq.3).or.(itimescheme.eq.7)) then ! AB3
               call decomp_2d_write_var(fh,disp,1,dphi1(:,:,:,3,is))
             end if
          end do
       endif
       call MPI_FILE_CLOSE(fh,ierror)
       ! Write info file for restart - Kay Schäfer
       if (nrank.eq.0) then
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
         write(111,'(A)')'/End'
         write(111,'(A)')'!========================='

         close(111)
       end if
    else
       if (nrank==0) then
         print *,'==========================================================='
         print *,'RESTART from file:', filestart
         print *,'==========================================================='
       end if
       call MPI_FILE_OPEN(MPI_COMM_WORLD, filestart, &
            MPI_MODE_RDONLY, MPI_INFO_NULL, &
            fh, ierror_o)
       disp = 0_MPI_OFFSET_KIND
       call decomp_2d_read_var(fh,disp,1,ux1)
       call decomp_2d_read_var(fh,disp,1,uy1)
       call decomp_2d_read_var(fh,disp,1,uz1)
       ! read previous time-step if necessary for AB2 or AB3
       if ((itimescheme.eq.2).or.(itimescheme.eq.3).or.(itimescheme.eq.7)) then ! AB2 or AB3
         call decomp_2d_read_var(fh,disp,1,dux1(:,:,:,2))
         call decomp_2d_read_var(fh,disp,1,duy1(:,:,:,2))
         call decomp_2d_read_var(fh,disp,1,duz1(:,:,:,2))
       end if
       ! for AB3 one more previous time-step
       if ((itimescheme.eq.3).or.(itimescheme.eq.7)) then ! AB3
         call decomp_2d_read_var(fh,disp,1,dux1(:,:,:,3))
         call decomp_2d_read_var(fh,disp,1,duy1(:,:,:,3))
         call decomp_2d_read_var(fh,disp,1,duz1(:,:,:,3))
       end if
       !
       call decomp_2d_read_var(fh,disp,3,pp3,phG)
       !
       if (iscalar==1) then
          do is=1, numscalar
             call decomp_2d_read_var(fh,disp,1,phi1(:,:,:,is))
             ! previous time-steps
             if ((itimescheme.eq.2).or.(itimescheme.eq.3).or.(itimescheme.eq.7)) then ! AB2 or AB3
               call decomp_2d_read_var(fh,disp,1,dphi1(:,:,:,2,is))
             end if
             !
             if ((itimescheme.eq.3).or.(itimescheme.eq.7)) then ! AB3
               call decomp_2d_read_var(fh,disp,1,dphi1(:,:,:,3,is))
             end if
          end do
       endif
       call MPI_FILE_CLOSE(fh,ierror_o)

       !! Read time of restart file
       write(filename,"('restart',I7.7,'.info')") ifirst-1
       inquire(file=filename, exist=fexists)
       if (nrank.eq.0) print *,filename
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

    if (nrank.eq.0) then
       if (ierror_o .ne. 0) then !Included by Felipe Schuch
          print *,'==========================================================='
          print *,'Error: Impossible to read '//trim(filestart)
          print *,'==========================================================='
          call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
       endif
    endif

    ! reconstruction of the dp/dx, dp/dy and dp/dz from pp3
    if (iresflg==0) then
       if (itimescheme.le.4) itr=1
       if (itimescheme.eq.5) itr=3
       if (itimescheme.eq.6) itr=5
       call gradp(px1,py1,pz1,pp3)
       if (nrank==0) print *,'reconstruction pressure gradients done!'
    end if

    if (iresflg .eq. 1 ) then !Writing restart
       if (nrank==0) then
          write(fmt1,"(I7.7)") itime
          print *,'Restart point restart',fmt1,' saved successfully!'!itime/icheckpoint,'saved successfully!'
          ! print *,'Elapsed time (s)',real(trestart,4)
          ! print *,'Aproximated writing speed (MB/s)',real(((s3df*16.)*1e-6)/trestart,4)
          print *,'If necesseary restart from:',itime+1
       endif
    end if

  end subroutine restart
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: channel_cfr
  !!      AUTHOR: Kay Schäfer
  !! DESCRIPTION: Inforces constant flow rate without need of data transposition
  !!
  !############################################################################
  subroutine channel_cfr (ux,constant)

    use decomp_2d
    use decomp_2d_poisson
    use variables
    use param
    use var
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux
    real(mytype) :: constant

    integer :: code,i,j,k,jloc
    real(mytype) :: can,ub,uball, dyloc
    !
    ub = zero
    uball = zero
    !
    do k=1,xsize(3)
       do j=xstart(2)+1,xend(2)-1
          jloc = j-xstart(2)+1
          dyloc  = (yp(j+1)-yp(j-1))
          do i=1,xsize(1)
            ub = ub + ux(i,jloc,k) * half * dyloc
          enddo
       enddo
    enddo

    ! Check if first and last index of subarray is at domain boundary
    if ( xstart(2)==1) then ! bottom point -> half distance
       ub = ub + sum(ux(:,1,:)) * yp(2)*half
    else
       ub = ub + sum(ux(:,1,:)) * (yp(xstart(2)+1)-yp(xstart(2)-1))*half
    end if
    !
    if (xend(2)==ny) then ! top point
       jloc = xend(2)-xstart(2)+1
       ub = ub + sum(ux(:,jloc,:)) * (yp(xend(2))-yp(xend(2)-1))*half
    else
       jloc = xend(2)-xstart(2)+1
       ub = ub + sum(ux(:,jloc,:)) * (yp(xend(2)+1)-yp(xend(2)-1))*half
    end if
    !
    ub = ub/(yly*(real(nx*nz,mytype)))

    call MPI_ALLREDUCE(ub,uball,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can=-(constant-uball)

    if (nrank==0) print *,nrank,'UT',uball,can

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux(i,j,k)=ux(i,j,k)-can
        enddo
      enddo
    enddo

    return
  end subroutine channel_cfr
  !############################################################################
  !##################################################################
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

     cfl_diff_x = xnu*dt/(dx**2)
     cfl_diff_z = xnu*dt/(dz**2)

     if (istret.eq.0) then
        cfl_diff_y   = xnu*dt/(dy**2)
     else
        cfl_diff_y = xnu*dt/(minval(dyp)**2)
     end if

     cfl_diff_sum = cfl_diff_x + cfl_diff_y + cfl_diff_z

     if (nrank==0) then
        print *,'==========================================================='
        print *,'Diffusion number'
        write(*,"(' cfl_diff_x             :        ',F13.8)") cfl_diff_x
        write(*,"(' cfl_diff_y             :        ',F13.8)") cfl_diff_y
        write(*,"(' cfl_diff_z             :        ',F13.8)") cfl_diff_z
        write(*,"(' cfl_diff_sum           :        ',F13.8)") cfl_diff_sum
        print *,'==========================================================='
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
    maxvalue_x  =-1609.
    maxvalue_y  =-1609.
    maxvalue_z  =-1609.
    maxvalue_sum=-1609.
    !
    if (istret.eq.0) then
       do j = xstart(2),xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:))/dx)
          value_y    = maxval(abs(uy(:,jloc,:))/dy)
          value_z    = maxval(abs(uz(:,jloc,:))/dz)
          value_sum  = maxval(abs(ux(:,jloc,:))/dx + abs(uy(:,jloc,:))/dy +    abs(uz(:,jloc,:))/dz)
          !
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    else
       do j = xstart(2),xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:))/dx)
          value_y    = maxval(abs(uy(:,jloc,:))/dyp(j))
          value_z    = maxval(abs(uz(:,jloc,:))/dz)
          value_sum  = maxval(abs(ux(:,jloc,:))/dx + abs(uy(:,jloc,:))/ dyp(j) + abs(uz(:,jloc,:))/dz)
          !
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    end if

    cflmax_in =  (/maxvalue_x, maxvalue_y, maxvalue_z, maxvalue_sum/)

    call    MPI_REDUCE(cflmax_in,cflmax_out,4,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    if (nrank.eq.0) then
      write(*,"(' CFL_x                  : ',F17.8)") cflmax_out(1)*dt
      write(*,"(' CFL_y                  : ',F17.8)") cflmax_out(2)*dt
      write(*,"(' CFL_z                  : ',F17.8)") cflmax_out(3)*dt
      !write(*,"(' CFL_sum                : ',F17.8)") cflmax_out(4)*dt
    end if
  end subroutine compute_cfl
  !##################################################################
  !##################################################################
  ! Rescale pressure to physical pressure
  ! Written by Kay Schäfer 2019
  !##################################################################
  subroutine rescale_pressure(pre1)

    use decomp_2d, only : nrank, mytype, xsize, ysize, zsize
    use param, only : itimescheme, gdt
    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: pre1

    ! Adjust pressure to physical pressure
    if  ((itimescheme.eq.2).or.(itimescheme.eq.3).or.(itimescheme.eq.5).or.(itimescheme.eq.7)) then !AB2, AB3, RK3, Semi-Impl. AB3
       pre1=pre1 / gdt(3) ! multiply pressure by factor of time-scheme (gdt = 1  / (dt * c_k) ) to get pyhsical pressure
    else
       if (nrank .eq. 0) print *,'WARNING: No scaling of pressure defined!!!'
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

    fm1 = sum(f1,DIM=1)/real(nx,mytype)
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

    fm2 = sum(f2,DIM=2)/real(ny,mytype)
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

    fm3 = sum(f3,DIM=3)/real(nz,mytype)
    return

  end subroutine mean_plane_z
end module tools
!##################################################################
!##################################################################
subroutine stabiltemp() !from Erik, adapted by Leonardo Romero Monteiro

  use param
  use variables
  use var

  implicit none

  complex(mytype) :: z,eit,ei2t,ei3t,eimt,eim2t,eim3t
  real(mytype) :: theta, dtheta, cc, fourier, cfl
  real(mytype) :: xkm, xk, xkp, xks, xkf, x, y
  real(mytype) :: am1, a0, a1, a2, a3
  real(mytype) :: bm1, b0, b1, b2, b3
  real(mytype) :: alpha1, c1, c11
  real(mytype) :: alpha2, c2
  real(mytype) :: alpha3, beta3, c3, d3
  integer :: i,ntheta,order

  ntheta=360
  dtheta=twopi/(ntheta-1.)
  xk=(fpi2+1.)*pi*pi
  order = 6   ! ordem da hiperviscosidade 0 = sem hiperviscosidade; 4 = 4a ordem com 2 formados;  6 = 6a ordem com 1 formado

  print *,'Writing stability data!'

  if (itimescheme==0) then !Euler (not implemented)
     am1=0; a0=1.; a1=0.; a2=0.
  endif

  if (itimescheme.eq.1) then !AB2
     am1=0; a0=1.5; a1=-0.5; a2=0.; a3=0.; bm1=1.; b0=-1.; b1=0.; b2=0.; b3=0.
  endif

  if (itimescheme.eq.3) then !RK3
     if (nrank==0) write(*,*) "Non implemented for RK3"
  endif

  if (itimescheme.eq.2) then !AB3
     am1=0.; a0=23./12.; a1=-16./12.; a2=5./12; a0=3./2+a2; a1=-1./2-2*a2; a3=0.; bm1=1.; b0=-1.; b1=0.; b2=0.; b3=0.
  endif

  open(10,file='stabiltemp_1.dat',form='formatted')
  do i=1,ntheta
     theta=(i-1)*dtheta

     eit=exp(cmplx(0.,1.)*theta)
     ei2t=eit*eit
     ei3t=eit*eit*eit
     eimt=1./eit
     eim2t=1./ei2t
     eim3t=1./ei3t
     !z=(eit-1.)/a0
     !z=(eit*(eit-1.))/(a0*eit+a1)
     !z=(ei3t-ei2t)/(a0*ei2t+a1*eit+a2)
     z=(bm1*eit+b0+b1*eimt+b2*eim2t+b3*eim3t)/(a0+a1*eimt+a2*eim2t+a3*eim3t)
     !z=(eit-1.)/(am1*eit+a0+a1*eimt)
     !z=(eit-1.)/(am1*eit+a0+a1*eimt+a2*eim2t)

     write(10,*) real(z),imag(z)
  enddo
  close(10)


  alpha1=1./3.
  a1=(alpha1+9.)/6.
  b1=(32.*alpha1-9.)/15.
  c1=(-3.*alpha1+1.)/10.

  if (order.eq.0) then

     alpha2=2./11
     a2=12./11
     b2=3./11
     c2=0.

  elseif (order.eq.4) then

     c11=exp(-((pi-2.*pi/3.)/(0.3*pi-2.*pi/3.))**2 )
     xkm=(c11*fpi2+1.)*(4./9.)*pi*pi

     alpha2=(64.*xkm-27.*xk-96.)/(64.*xkm-54.*xk+48.)
     a2 = (54.*xk-15.*xkm*xk+12.)/(64.*xkm-54.*xk+48.)
     b2 = (192.*xkm-216.*xk+24.*xkm*xk-48.)/(64.*xkm-54.*xk+48.)
     c2 = 3.*(18.*xk -3.*xkm*xk-36.)/(64.*xkm-54.*xk+48.)

  elseif(order.eq.6) then

     alpha2=(45.*xk-272.)/(2*(45.*xk-208.))
     c2=(2.-11.*alpha2)/20.
     a2=(6.-9.*alpha2)/4.
     b2=(-3.+24.*alpha2)/5.

  endif

  !alpha3=0.45
  !beta3=(3.-2.*alpha3)/10.
  !a3=(2.+3.*alpha3)/4.
  !b3=(6.+7*alpha3)/8.
  !c3=(6.+alpha3)/20.
  !d3=(2-3.*alpha3)/40.

  cc=4.
  fourier=xnu*dt/(dx*dx)
  cfl=cc*dt/dx

  open(10,file='stabiltemp_2.dat',form='formatted')
  do i=1,ntheta
     theta=(i-1)*dtheta

     xkp=(a1*sin(theta)+(b1/2)*sin(2*theta) +(c1/3)*sin(3*theta))/(1+2*alpha1*cos(theta))
     xks=(2*a2*(1-cos(theta))+(b2/2)*(1-cos(2*theta)) +(2*c2/9)*(1-cos(3*theta)))/(1+2*alpha2*cos(theta))
     !xkf=(a3+b3*cos(theta)+c3*cos(2*theta)+d3*cos(3*theta)) /(1+2*alpha3*cos(theta)+2*beta3*cos(2*theta))
     x=-fourier*xks
     y=-cfl*xkp!*xkf

     write(10,*) x,y
  enddo
  close(10)

end subroutine stabiltemp
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

  implicit none

  real(mytype),intent(in) :: uxmax,uymax,uzmax
  real(mytype) :: cfl_x_adv,cfl_x_diff,cfl_y_adv,cfl_y_diff,cfl_z_adv,cfl_z_diff
  real(mytype) :: cfl_conv_lim, cfl_diff_lim
  real(mytype) :: sigma_conv(3), sigma_diff(3)
  real(mytype) :: visc

  ! Set the constants (this is true for periodic boundaries)
  sigma_conv=[0.0, sqrt(3.0), 2.85]
  sigma_diff=[2.0, 2.5, 2.9]

  if(jles==0) then
     visc=xnu
  elseif (jles==1) then
     visc=20*fpi2*xnu
  endif

  ! This is considering 1D peridic boundaries
  ! Do x-direction
  cfl_x_adv=abs(uxmax)*dt/dx
  cfl_x_diff=visc*dt/dx**2.0
  ! Do y-direction
  cfl_y_adv=abs(uymax)*dt/dy
  cfl_y_diff=visc*dt/dy**2.0
  ! Do z-direction
  cfl_z_adv=abs(uzmax)*dt/dz
  cfl_z_diff=visc*dt/dz**2.0

  ! So far we will focus on uniform grids
  if(nrank==0) then
     write(*,*) ' '
     write(*,1002) cfl_x_adv, cfl_x_diff
1002 format('CFL x-direction (Adv and Diff) =',F9.4,',',F9.4)
     write(*,1003) cfl_y_adv, cfl_y_diff
1003 format('CFL y-direction (Adv and Diff) =',F9.4,',',F9.4)
     write(*,1004) cfl_z_adv, cfl_z_diff
1004 format('CFL z-direction (Adv and Diff) =',F9.4,',',F9.4)
     cfl_conv_lim=sigma_conv(itimescheme)/sqrt(3.0)
     cfl_diff_lim=sigma_diff(itimescheme)/6.0
     write(*,1005) cfl_conv_lim, cfl_diff_lim
     write(*,*) ' '
1005 format('CFL limits (Adv and Diff) : ',F9.4,',',F9.4)
  endif

end subroutine cfl_compute
!##################################################################
!##################################################################
subroutine stretching()

  USE decomp_2d
  !USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

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

  if (nrank==0) then
     open(10,file='yp.dat', form='formatted')
     do j=1,ny
        write(10,*)yp(j)
     enddo
     close(10)
  endif

end subroutine stretching
!##################################################################
!##################################################################
subroutine inversion5_v1(aaa,eee,spI)

  USE decomp_2d
  !USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  do i=1,2
     ja(i)=4-i
     jb(i)=5-i
  enddo
  do m=1,ny/2-2
     do i=1,2
        mi=m+i
        do k=spI%yst(3),spI%yen(3)
           do j=spI%yst(1),spI%yen(1)
              if (real(aaa(j,m,k,3), kind=mytype).ne.zero) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
              if (aimag(aaa(j,m,k,3)).ne.zero)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
              sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
              eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
                   aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
           enddo
        enddo
        do jc=ja(i),jb(i)
           do k=spI%yst(3),spI%yen(3)
              do j=spI%yst(1),spI%yen(1)
                 aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
                      aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
              enddo
           enddo
        enddo
     enddo
  enddo


  do k=spI%yst(3),spI%yen(3)
     do j=spI%yst(1),spI%yen(1)
        if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
           tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
        else
           tmp1=zero
        endif
        if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
           tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
        else
           tmp2=zero
        endif
        sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
             aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

        if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
           tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
           tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), kind=mytype)
        else
           tmp1=zero
           tmp3=zero
        endif
        if (abs(aimag(b1(j,k))).gt.epsilon) then
           tmp2=aimag(sr(j,k))/aimag(b1(j,k))
           tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
        else
           tmp2=zero
           tmp4=zero
        endif
        a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)

        if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
           tmp1=one/real(aaa(j,ny/2-1,k,3), kind=mytype)
        else
           tmp1=zero
        endif
        if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
           tmp2=one/aimag(aaa(j,ny/2-1,k,3))
        else
           tmp2=zero
        endif
        b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
             aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
        eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
             aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
     enddo
  enddo

  do i=ny/2-2,1,-1
     do k=spI%yst(3),spI%yen(3)
        do j=spI%yst(1),spI%yen(1)
           if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
              tmp1=one/real(aaa(j,i,k,3), kind=mytype)
           else
              tmp1=zero
           endif
           if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
              tmp2=one/aimag(aaa(j,i,k,3))
           else
              tmp2=zero
           endif
           sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
           a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
                aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
           b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
                aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
           eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
                real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
                real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
                aimag(eee(j,i,k))*aimag(sr(j,k))-&
                aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v1
!##################################################################
!##################################################################
subroutine inversion5_v2(aaa,eee,spI)

  USE decomp_2d
  !USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

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

  do i=1,2
     ja(i)=4-i
     jb(i)=5-i
  enddo
  do m=1,nym-2
     do i=1,2
        mi=m+i
        do k=spI%yst(3),spI%yen(3)
           do j=spI%yst(1),spI%yen(1)
              if (real(aaa(j,m,k,3), kind=mytype).ne.zero) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
              if (aimag(aaa(j,m,k,3)).ne.zero)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
              sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
              eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
                   aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
           enddo
        enddo
        do jc=ja(i),jb(i)
           do k=spI%yst(3),spI%yen(3)
              do j=spI%yst(1),spI%yen(1)
                 aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
                      aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
              enddo
           enddo
        enddo
     enddo
  enddo
  do k=spI%yst(3),spI%yen(3)
     do j=spI%yst(1),spI%yen(1)
        if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
           tmp1=real(aaa(j,nym,k,2), kind=mytype)/real(aaa(j,nym-1,k,3), kind=mytype)
        else
           tmp1=zero
        endif
        if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
           tmp2=aimag(aaa(j,nym,k,2))/aimag(aaa(j,nym-1,k,3))
        else
           tmp2=zero
        endif
        sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        b1(j,k)=cmplx(real(aaa(j,nym,k,3), kind=mytype)-tmp1*real(aaa(j,nym-1,k,4), kind=mytype),&
             aimag(aaa(j,nym,k,3))-tmp2*aimag(aaa(j,nym-1,k,4)), kind=mytype)
        if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
           tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
           tmp3=real(eee(j,nym,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,nym-1,k), kind=mytype)
        else
           tmp1=zero
           tmp3=zero
        endif
        if (abs(aimag(b1(j,k))).gt.epsilon) then
           tmp2=aimag(sr(j,k))/aimag(b1(j,k))
           tmp4=aimag(eee(j,nym,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,nym-1,k))
        else
           tmp2=zero
           tmp4=zero
        endif
        a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        eee(j,nym,k)=cmplx(tmp3,tmp4, kind=mytype)

        if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
           tmp1=one/real(aaa(j,nym-1,k,3), kind=mytype)
        else
           tmp1=zero
        endif
        if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
           tmp2=one/aimag(aaa(j,nym-1,k,3))
        else
           tmp2=zero
        endif
        b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        a1(j,k)=cmplx(real(aaa(j,nym-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
             aimag(aaa(j,nym-1,k,4))*aimag(b1(j,k)), kind=mytype)
        eee(j,nym-1,k)=cmplx(real(eee(j,nym-1,k), kind=mytype)*real(b1(j,k), kind=mytype)-&
             real(a1(j,k), kind=mytype)*real(eee(j,nym,k), kind=mytype),&
             aimag(eee(j,nym-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,nym,k)), kind=mytype)
     enddo
  enddo

  do i=nym-2,1,-1
     do k=spI%yst(3),spI%yen(3)
        do j=spI%yst(1),spI%yen(1)
           if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
              tmp1=one/real(aaa(j,i,k,3), kind=mytype)
           else
              tmp1=zero
           endif
           if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
              tmp2=one/aimag(aaa(j,i,k,3))
           else
              tmp2=zero
           endif
           sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
           a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
                aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
           b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
                aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
           eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
                real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
                real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
                aimag(eee(j,i,k))*aimag(sr(j,k))-&
                aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v2
!##################################################################
!##################################################################
subroutine tripping(tb,ta)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

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

  if ((itime.eq.ifirst).and.(nrank.eq.0)) then
     call random_seed(SIZE=ii)
     call random_seed(PUT=seed0*(/ (1, i = 1, ii) /))

     !DEBUG:
     !call random_number(randx)
     !call MPI_BCAST(randx,1,real_type,0,MPI_COMM_WORLD,code)
     !write(*,*) 'RANDOM:', nrank, randx, ii
     !First random generation of h_nxt


     do j=1,z_modes

        call random_number(randx)
        h_coeff(j)=1.0*(randx-0.5)
     enddo
     h_coeff=h_coeff/sqrt(DBLE(z_modes))
  endif

  !Initialization h_nxt  (always bounded by xsize(3)^2 operations)
  if (itime.eq.ifirst) then
     call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
     nxt_itr=0
     do k=1,xsize(3)
        h_nxt(k)=0.0
        z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(2.0*pi*j*z_pos/zlz)
        enddo
     enddo
  end if



  !Time-loop
  i=int(t/ts_tr)
  if (i.ge.nxt_itr) then  !Nxt_itr is a global variable
     nxt_itr=i+1

     !First random generation of h
     h_i(:)=h_nxt(:)
     if (nrank .eq. 0) then
        do j=1,z_modes
           call random_number(randx)
           h_coeff(j)=1.0*(randx-0.5)
        enddo
        h_coeff=h_coeff/sqrt(DBLE(z_modes)) !Non-dimensionalization
     end if

     call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)


     !Initialization h_nxt  (always bounded by z_steps^2 operations)
     do k=1,xsize(3)
        h_nxt(k)=0.0
        z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(2.0*pi*j*z_pos/zlz)
        enddo
     enddo
  endif

  !Time coefficient
  p_tr=t/ts_tr-i
  b_tr=3.0*p_tr**2-2.0*p_tr**3

  !Creation of tripping velocity
  do i=1,xsize(1)
     x_pos=(xstart(1)+(i-1)-1)*dx
     do j=1,xsize(2)
        !y_pos=(xstart(2)+(j-1)-1)*dy
        y_pos=yp(xstart(2)+(j-1))
        do k=1,xsize(3)
           !g(z)*EXP_F(X,Y)
           ta(i,j,k)=((1.0-b_tr)*h_i(k)+b_tr*h_nxt(k))
           !ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-(y_pos/ys_tr)**2)*ta(i,j,k)
           ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-((y_pos-0.5)/ys_tr)**2)*ta(i,j,k)
           tb(i,j,k)=tb(i,j,k)+ta(i,j,k)

           z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
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

  USE param
  USE variables
  USE decomp_2d
  USE MPI

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


  if ((itime.eq.ifirst).and.(nrank.eq.0)) then
     call random_seed(SIZE=ii)
     call random_seed(PUT=seed0*(/ (1, i = 1, ii) /))

     INQUIRE(FILE='restart.nc',exist=exist)
     !if ((ilit==1).AND.(exist)) then
     !if (exist) then
     !print*, 'h_coeff1 and phase1 already read from restart.nc'
     !print*, 'h_coeff2 and phase2 already read from restart.nc'
     !nxt_itr=int(t/ts_tr_tbl)
     !else
     nxt_itr=1
     do j=1,z_modes
        call random_number(randx)
        h_coeff1(j)=1.0*(randx-0.5)/sqrt(DBLE(z_modes))
        call random_number(randx)
        phase1(j) = 2.0*pi*randx
        call random_number(randx)
        h_coeff2(j)=1.0*(randx-0.5)/sqrt(DBLE(z_modes))
        call random_number(randx)
        phase2(j) = 2.0*pi*randx
     enddo
     !endif
  endif

  !Initialization h_nxt  (always bounded by xsize(3)^2 operations)
  if (itime.eq.ifirst) then
     call MPI_BCAST(h_coeff1,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(phase1,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(h_coeff2,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(phase2,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(nxt_itr,1,mpi_int,0,MPI_COMM_WORLD,code)

     do k=1,xsize(3)
        h_1(k)=0.0
        h_2(k)=0.0
        z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_1(k)= h_1(k)+h_coeff1(j)*sin(2.0*pi*j*z_pos/zlz+phase1(j))
           h_2(k)= h_2(k)+h_coeff2(j)*sin(2.0*pi*j*z_pos/zlz+phase2(j))
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
     if (nrank .eq. 0) then
        do j=1,z_modes
           call random_number(randx)
           h_coeff1(j)=1.0*(randx-0.5)/sqrt(DBLE(z_modes))
           call random_number(randx)
           phase1(j) = 2.0*pi*randx
        enddo
     end if

     call MPI_BCAST(h_coeff1,z_modes,real_type,0,MPI_COMM_WORLD,code)
     call MPI_BCAST(phase1,z_modes,real_type,0,MPI_COMM_WORLD,code)

     !Initialization h_nxt  (always bounded by z_steps^2 operations)
     do k=1,xsize(3)
        h_1(k)=0.0
        z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
        do j=1,z_modes
           h_1(k)= h_1(k)+h_coeff1(j)*sin(2.0*pi*j*z_pos/zlz+phase1(j))
        enddo
     enddo
  endif
  !-------------------------------------------------------------------

  !Time coefficient
  p_tr=t/ts_tr_tbl-i
  b_tr=3.0*p_tr**2-2.0*p_tr**3
  !Creation of tripping velocity
  do i=1,xsize(1)
     x_pos=(xstart(1)+(i-1)-1)*dx
     do j=1,xsize(2)
        y_pos=yp(xstart(2)+(j-1))
        do k=1,xsize(3)
           ta(i,j,k)=((1.0-b_tr)*h_1(k)+b_tr*h_2(k))
           ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr_tbl)/xs_tr_tbl)**2-((y_pos-0.05)/ys_tr_tbl)**2)*ta(i,j,k)
           tb(i,j,k)=tb(i,j,k)+ta(i,j,k)

           z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz

        enddo
     enddo
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,code)
  !if (nrank==0) print*, maxval(ta(:,:,:)),minval(ta), z_modes

  return
end subroutine tbl_tripping
!##################################################################
!##################################################################
function rl(complexnumber)

  USE param

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!##################################################################
!##################################################################
function iy(complexnumber)

  USE param

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!##################################################################
!##################################################################
function cx(realpart,imaginarypart)

  USE param

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!##################################################################
!##################################################################
SUBROUTINE calc_temp_eos(temp, rho, phi, mweight, xlen, ylen, zlen)

  USE decomp_2d
  USE param, ONLY : pressure0, imultispecies
  USE var, ONLY : numscalar

  IMPLICIT NONE

  !! INPUTS
  INTEGER, INTENT(IN) :: xlen, ylen, zlen
  REAL(mytype), INTENT(IN), DIMENSION(xlen, ylen, zlen) :: rho
  REAL(mytype), INTENT(IN), DIMENSION(xlen, ylen, zlen, numscalar) :: phi

  !! OUTPUTS
  REAL(mytype), INTENT(OUT), DIMENSION(xlen, ylen, zlen) :: temp

  !! LOCALS
  REAL(mytype), DIMENSION(xlen, ylen, zlen) :: mweight

  temp(:,:,:) = pressure0 / rho(:,:,:)
  IF (imultispecies) THEN
     CALL calc_mweight(mweight, phi, xlen, ylen, zlen)
     temp(:,:,:) = temp(:,:,:) * mweight(:,:,:)
  ENDIF

ENDSUBROUTINE calc_temp_eos
!##################################################################
!##################################################################
SUBROUTINE calc_rho_eos(rho, temp, phi, mweight, xlen, ylen, zlen)

  USE decomp_2d
  USE param, ONLY : pressure0, imultispecies
  USE var, ONLY : numscalar

  IMPLICIT NONE

  !! INPUTS
  INTEGER, INTENT(IN) :: xlen, ylen, zlen
  REAL(mytype), INTENT(IN), DIMENSION(xlen, ylen, zlen) :: temp
  REAL(mytype), INTENT(IN), DIMENSION(xlen, ylen, zlen, numscalar) :: phi

  !! OUTPUTS
  REAL(mytype), INTENT(OUT), DIMENSION(xlen, ylen, zlen) :: rho

  !! LOCALS
  REAL(mytype), DIMENSION(xlen, ylen, zlen) :: mweight

  rho(:,:,:) = pressure0 / temp(:,:,:)
  IF (imultispecies) THEN
     CALL calc_mweight(mweight, phi, xlen, ylen, zlen)
     rho(:,:,:) = rho(:,:,:) * mweight(:,:,:)
  ENDIF

ENDSUBROUTINE calc_rho_eos
!##################################################################
!##################################################################
SUBROUTINE calc_mweight(mweight, phi, xlen, ylen, zlen)

  USE decomp_2d
  USE param, ONLY : zero, one
  USE param, ONLY : massfrac, mol_weight
  USE var, ONLY : numscalar

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: xlen, ylen, zlen
  REAL(mytype), INTENT(IN), DIMENSION(xlen, ylen, zlen, numscalar) :: phi

  !! LOCALS
  REAL(mytype), DIMENSION(xlen, ylen, zlen) :: mweight
  INTEGER :: is

  mweight(:,:,:) = zero
  DO is = 1, numscalar
     IF (massfrac(is)) THEN
        mweight(:,:,:) = mweight(:,:,:) + phi(:,:,:,is) / mol_weight(is)
     ENDIF
  ENDDO
  mweight(:,:,:) = one / mweight(:,:,:)

ENDSUBROUTINE calc_mweight
!##################################################################
!##################################################################
