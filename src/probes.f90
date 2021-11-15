!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: probes.f90
!!!      AUTHOR: Felipe Schuch
!!!    MODIFIED: Cédric Flageul
!!! DESCRIPTION: This module is dedicated to monitoring points.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module probes

  USE decomp_2d, only : ph1, nrank, mytype
  USE decomp_2d, only : xstart, xend, ystart, yend, zstart, zend
  USE decomp_2d, only : decomp_2d_abort

  IMPLICIT NONE

  !
  ! Following parameters are extracted from the i3d input file
  !
  ! Number of probes
  integer, save :: nprobes
  ! Flag to monitor velocity, scalar(s) and pressure gradients
  logical, save :: flag_extra_probes = .false.
  ! Flag to print 16 digits, only 6 digits by default
  logical, save :: flag_all_digits = .false.
  ! Position of the probes
  real(mytype), save, allocatable, dimension(:,:) :: xyzprobes
  !
  ! Following parameters are not extracted from the i3d input file
  !
  ! Offset in case of restart without probes
  integer, save :: main_probes_offset, extra_probes_offset
  ! Probes location on the velocity grid
  logical, save, allocatable, dimension(:) :: rankprobes, rankprobesY, rankprobesZ
  integer, save, allocatable, dimension(:) :: nxprobes, nyprobes, nzprobes
  ! Probes location on the pressure grid
  logical, save, allocatable, dimension(:) :: rankprobesP
  integer, save, allocatable, dimension(:) :: nxprobesP, nyprobesP, nzprobesP
  ! Arrays used to monitor gradients
  real(mytype), save, allocatable, dimension(:,:) :: dfdx, dfdy, dfdz

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: setup_probes, init_probes, write_probes, finalize_probes, &
            nprobes, flag_all_digits, flag_extra_probes, xyzprobes

contains

  !############################################################################
  !
  ! This subroutine is used to setup the probes module, before init_probes
  !
  subroutine setup_probes()

    ! Exit if no file or no probes
    if (nprobes.le.0) then
       flag_extra_probes = .false.
       return
    endif

    ! Allocate memory
    allocate(xyzprobes(3, nprobes))

  end subroutine

  !############################################################################
  !
  ! This subroutine is used to initialize the probes module, after setup_probes
  !
  subroutine init_probes()

    USE decomp_2d, only : real_type
    USE MPI
    USE param, only : dx, dy, dz, nclx, ncly, nclz, xlx, yly, zlz, istret, one, half
    USE param, only : irestart, ifirst
    USE variables, only : nxm, ny, yp, nym, ypi, nzm, numscalar

    logical :: fexists
    integer :: i,j,code
    real(mytype) :: xprobes, yprobes, zprobes
    character(len=30) :: filename

    if (nprobes.le.0) return

#ifdef DEBG
    if (nrank == 0) print *,'# init_probes start'
#endif

    ! In case of restart, check existence of previous probes results
    main_probes_offset = 0
    extra_probes_offset = 0
    !
    if (irestart.ne.0) then

       ! Basic probes
       if (nrank.eq.0) then
          write(filename,"('./probes/probe',I4.4)") 1
          inquire(file=filename, exist=fexists)
          if (.not.fexists) then
             main_probes_offset = ifirst - 1
          endif
       endif
       ! Broadcast
       call MPI_BCAST(main_probes_offset,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)

       ! Extra probes
       if (flag_extra_probes) then
          if (nrank.eq.0) then
             write(filename,"('./probes/probe_dx_',I4.4)") 1
             inquire(file=filename, exist=fexists)
             if (.not.fexists) then
                extra_probes_offset = ifirst - 1
             endif
          endif
          ! Broadcast
          call MPI_BCAST(extra_probes_offset,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
       endif

    endif

    ! Probes on the velocity grid
    allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes))
    allocate(rankprobes(nprobes), rankprobesY(nprobes), rankprobesZ(nprobes))
    rankprobes(:) = .false.
    rankprobesY(:) = .false.
    rankprobesZ(:) = .false.
    ! Probes on the pressure grid
    allocate(nxprobesP(nprobes), nyprobesP(nprobes), nzprobesP(nprobes))
    allocate(rankprobesP(nprobes))
    rankprobesP(:) = .false.

    do i = 1, nprobes

       ! Enforce bounds
       xprobes = max(epsilon(xlx), min(xlx*(one-epsilon(xlx)), xyzprobes(1, i)))
       yprobes = max(epsilon(xlx), min(yly*(one-epsilon(xlx)), xyzprobes(2, i)))
       zprobes = max(epsilon(xlx), min(zlz*(one-epsilon(xlx)), xyzprobes(3, i)))

       !x
       ! 1 <= nxprobes(i) <= nx
       ! x = (nxprobes(i) - 1) * Lx / (nx-1)
       if (nclx) then
          nxprobes(i) = int(xprobes/dx) + 1
       else
          nxprobes(i) = nint(xprobes/dx) + 1
       endif
       ! 1 <= nxprobesP(i) <= nxm
       ! xp = (nxprobesP(i) - 1 + 0.5) * Lx / nxm
       nxprobesP(i) = nint(nxm*xprobes/xlx-half) + 1

       !y
       if (istret.eq.0) then
          if (ncly) then
             nyprobes(i) = int(yprobes/dy) + 1
          else
             nyprobes(i) = nint(yprobes/dy) + 1
          endif
          nyprobesP(i) = nint(nym*yprobes/yly-half) + 1
       else
          if (yprobes.le.ypi(1)) then
             nyprobes(i) = 1
          else if (yprobes.ge.ypi(nym)) then
             nyprobes(i) = ny
          else
             do j = 1, nym-1
                if (ypi(j).le.yprobes .and. yprobes.lt.ypi(j+1)) then
                   nyprobes(i) = j + 1
                end if
             end do
          endif
          if (yprobes.le.yp(2)) then
             nyprobesP(i) = 1
          elseif (yprobes.ge.yp(ny-1)) then
             nyprobesP(i) = nym
          else
             do j = 3, ny-2
                if (yp(j).le.yprobes .and. yprobes.lt.yp(j+1)) then
                   nyprobesP(i) = j - 1
                   exit
                end if
             end do
          end if
       end if

       !z
       if (nclz) then
          nzprobes(i) = int(zprobes/dz) + 1
       else
          nzprobes(i) = nint(zprobes/dz) + 1
       endif
       nzprobesP(i) = nint(nzm*zprobes/zlz-half) + 1

       ! Flag the rank with the probe
       if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
          if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
             if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                rankprobes(i) = .true.
             endif
          endif
       endif
       if       (ystart(1) .le. nxprobes(i) .and. nxprobes(i) .le. yend(1)) then
          if    (ystart(2) .le. nyprobes(i) .and. nyprobes(i) .le. yend(2)) then
             if (ystart(3) .le. nzprobes(i) .and. nzprobes(i) .le. yend(3)) then
                rankprobesY(i) = .true.
             endif
          endif
       endif
       if       (zstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. zend(1)) then
          if    (zstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. zend(2)) then
             if (zstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. zend(3)) then
                rankprobesZ(i) = .true.
             endif
          endif
       endif
       if       (ph1%zst(1) .le. nxprobesP(i) .and. nxprobesP(i) .le. ph1%zen(1)) then
          if    (ph1%zst(2) .le. nyprobesP(i) .and. nyprobesP(i) .le. ph1%zen(2)) then
             if (ph1%zst(3) .le. nzprobesP(i) .and. nzprobesP(i) .le. ph1%zen(3)) then
                rankprobesP(i) = .true.
             endif
          endif
       endif
    enddo

    ! Log the real position of the probes in a file
    if (nrank.eq.0) call log_probes_position()

    ! Allocate memory to monitor gradients
    if (flag_extra_probes) then
      ! X gradients : velocity + scalar(s) + 3 x pressure
      allocate(dfdx(nprobes, 3+numscalar+3))
      ! Y gradients : velocity + scalar(s)
      allocate(dfdy(nprobes, 3+numscalar))
      ! Z gradients : velocity + scalar(s)
      allocate(dfdz(nprobes, 3+numscalar))
    endif

#ifdef DEBG
    if (nrank == 0) print *,'# init_probes ok'
#endif

  end subroutine init_probes

  !############################################################################
  !
  ! This subroutine is used to monitor velocity, scalar(s) and pressure
  !
  subroutine write_probes(ux1,uy1,uz1,pp3,phi1)

    use decomp_2d, only : xsize, ysize, zsize
    use param, only : itime, irestart, itr, t
    use param, only : npress, sync_vel_needed, sync_scal_needed
    use var, only : nzmsize
    use variables, only : numscalar
    use variables, only : derx, dery, derz, derxs, derys, derzs
    use var, only : transpose_x_to_y, transpose_y_to_z
    use var, only : di1, di2, di3
    use var, only : ffx, ffxp, ffxS, ffxpS, fwxpS, fsx, fsxp, fsxpS, fsxS, fwx, fwxp, fwxS, sx
    use var, only : ffy, ffyp, ffyS, ffypS, fwypS, fsy, fsyp, fsypS, fsyS, fwy, fwyp, fwyS, sy, ppy
    use var, only : ffz, ffzp, ffzS, ffzpS, fwzpS, fsz, fszp, fszpS, fszS, fwz, fwzp, fwzS, sz

    use var, only : zero
    use var, only : ux2, uy2, uz2, phi2
    use var, only : ux3, uy3, uz3, phi3
    use var, only : px1, py1, pz1
    use var, only : ta1, tb1, tc1, td1, te1, tf1
    use var, only : tc2, td2, te2, tf2
    use var, only : tb3, td3, te3, tf3
    use var, only : sc_even
    
    use ibm_param, only : ubcx,ubcy,ubcz

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype), intent(in),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

    integer :: iounit, i, FS, FSP, is
    character(len=100) :: fileformat, fileformatP
    character(len=1),parameter :: NL=char(10) !new line character
    character(len=30) :: filename
    logical :: evensc
    
    if (nprobes.le.0) return

    ! Number of columns
    FS = 1+3+numscalar
    FSP = 1+1 ! Pressure grid
    ! All digits or only 6 digits
    if (flag_all_digits) then
       write(fileformat, '( "(",I4,"(E24.16),A)" )' ) FS
       write(fileformatP, '( "(",I4,"(E24.16),A)" )' ) FSP
       ! Line width
       FS = FS*24+1
       FSP = FSP*24+1
    else
       write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
       write(fileformatP, '( "(",I4,"(E14.6),A)" )' ) FSP
       ! Line width
       FS = FS*14+1
       FSP = FSP*14+1
    endif

    do i=1, nprobes
       if (rankprobes(i)) then
          write(filename,"('./probes/probe',I4.4)") i
          open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          if (numscalar.ge.1) then
             write(iounit,fileformat,rec=itime-main_probes_offset) t,&  !1
                  ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
                  uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
                  uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
                  phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !numscalar
                  NL                                                    !+1
          else
             write(iounit,fileformat,rec=itime-main_probes_offset) t,&  !1
                  ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
                  uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
                  uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
                  NL                                                    !+1
          endif
          close(iounit)
       endif
       if (rankprobesP(i)) then
          write(filename,"('./probes/probeP',I4.4)") i
          open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FSP)
          write(iounit,fileformatP,rec=itime-main_probes_offset) t,&    !1
               pp3(nxprobesP(i),nyprobesP(i),nzprobesP(i),1),&          !2
               NL                                                       !+1
          close(iounit)
       endif
    enddo

    ! Monitor gradients
    if (flag_extra_probes) then

       ! Perform comunications if needed
       if (sync_vel_needed) then
          call transpose_x_to_y(ux1,ux2)
          call transpose_x_to_y(uy1,uy2)
          call transpose_x_to_y(uz1,uz2)
          call transpose_y_to_z(ux2,ux3)
          call transpose_y_to_z(uy2,uy3)
          call transpose_y_to_z(uz2,uz3)
          sync_vel_needed = .false.
       endif
       ! Compute velocity gradient
       call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
       call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
       call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
       call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
       call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
       call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
       call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
       call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
       call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
       ! Store velocity gradient
       call write_extra_probes_vel(ta1, td2, td3, tb1, te2, te3, tc1, tf2, tf3)

       ! Store pressure gradient
       call write_extra_probes_pre(px1, py1, pz1)

       ! Monitor scalars gradient 
       do is = 1, numscalar
          evensc = .true.
          if (.not.sc_even(is)) then
             evensc = .false.
          endif

          ! Perform communications if needed
          if (sync_scal_needed) then
             call transpose_x_to_y(phi1(:,:,:,is),phi2(:,:,:,is))
             call transpose_y_to_z(phi2(:,:,:,is),phi3(:,:,:,is))
          endif
          ! Compute derivative
          if (evensc) then
             call derxS (tb1,phi1(:,:,:,is),di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1,zero)
             call deryS (tc2,phi2(:,:,:,is),di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1,zero)
             call derzS (tb3,phi3(:,:,:,is),di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1,zero)
          else
             call derxS (tb1,phi1(:,:,:,is),di1,sx,ffxS,fsxS,fwxS,xsize(1),xsize(2),xsize(3),0,zero)
             call deryS (tc2,phi2(:,:,:,is),di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),0,zero)
             call derzS (tb3,phi3(:,:,:,is),di3,sz,ffzS,fszS,fwzS,zsize(1),zsize(2),zsize(3),0,zero)
          endif
          ! Store scalars gradient
          call write_extra_probes_scal(is, tb1, tc2, tb3)

       end do

       ! Scalars are synchronized
       sync_scal_needed = .false.

       ! Record the gradients
       call write_extra_probes()

    endif

  end subroutine write_probes

  !############################################################################
  !
  ! This subroutine is used to monitor velocity gradient
  !
  subroutine write_extra_probes_vel(duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz)

    ! Arguments
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: duxdx, duydx, duzdx
    real(mytype),intent(in),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: duxdy, duydy, duzdy
    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: duxdz, duydz, duzdz

    ! Local variables
    integer :: i

    do i = 1, nprobes
       if (rankprobes(i)) then
          dfdx(i,1) = duxdx(nxprobes(i),nyprobes(i),nzprobes(i))
          dfdx(i,2) = duydx(nxprobes(i),nyprobes(i),nzprobes(i))
          dfdx(i,3) = duzdx(nxprobes(i),nyprobes(i),nzprobes(i))
       endif
       if (rankprobesY(i)) then
          dfdy(i,1) = duxdy(nxprobes(i),nyprobes(i),nzprobes(i))
          dfdy(i,2) = duydy(nxprobes(i),nyprobes(i),nzprobes(i))
          dfdy(i,3) = duzdy(nxprobes(i),nyprobes(i),nzprobes(i))
       endif
       if (rankprobesZ(i)) then
          dfdz(i,1) = duxdz(nxprobes(i),nyprobes(i),nzprobes(i))
          dfdz(i,2) = duydz(nxprobes(i),nyprobes(i),nzprobes(i))
          dfdz(i,3) = duzdz(nxprobes(i),nyprobes(i),nzprobes(i))
       endif
    enddo

  end subroutine write_extra_probes_vel

  !############################################################################
  !
  ! This subroutine is used to monitor scalar(s) gradient
  !
  subroutine write_extra_probes_scal(is, dphidx, dphidy, dphidz)

    use variables, only : numscalar

    ! Arguments
    integer, intent(in) :: is
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dphidx
    real(mytype),intent(in),dimension(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)) :: dphidy
    real(mytype),intent(in),dimension(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)) :: dphidz

    ! Local variables
    integer :: i

    do i = 1, nprobes
       if (rankprobes(i)) then
          dfdx(i,3+is) = dphidx(nxprobes(i),nyprobes(i),nzprobes(i))
       endif
       if (rankprobesY(i)) then                                                                       
          dfdy(i,3+is) = dphidy(nxprobes(i),nyprobes(i),nzprobes(i))
       endif
       if (rankprobesZ(i)) then                                                                       
          dfdz(i,3+is) = dphidz(nxprobes(i),nyprobes(i),nzprobes(i))
       endif
    enddo

  end subroutine write_extra_probes_scal

  !############################################################################
  !
  ! This subroutine is used to monitor the pressure gradient
  !
  subroutine write_extra_probes_pre(dpdx, dpdy, dpdz)

    use param, only : one, gdt, itimescheme
    use variables, only : numscalar

    ! Arguments
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: dpdx, dpdy, dpdz

    ! Local variables
    integer :: i
    real(mytype) :: fact

    ! Explicit Euler, AB2, AB3, AB4, RK3
    if (itimescheme.ge.1 .and. itimescheme.le.5) then
       fact = one / gdt(3)
    ! RK4
    elseif (itimescheme.eq.6) then
       fact = one / gdt(5)
    else
       ! We should not be here
       fact = one
    endif

    do i = 1, nprobes
       if (rankprobes(i)) then
          dfdx(i,3+numscalar+1) = dpdx(nxprobes(i),nyprobes(i),nzprobes(i)) * fact
          dfdx(i,3+numscalar+2) = dpdy(nxprobes(i),nyprobes(i),nzprobes(i)) * fact
          dfdx(i,3+numscalar+3) = dpdz(nxprobes(i),nyprobes(i),nzprobes(i)) * fact
       endif
    enddo

  end subroutine write_extra_probes_pre

  !############################################################################
  !
  ! This subroutine performs IO.
  ! Gradients are recorded at the beginning of each time step when itr = 1.
  !
  subroutine write_extra_probes()

    use param, only : itime, t, dt
    use variables, only : numscalar

    integer :: iounit, i, FSX, FSY, FSZ
    character(len=100) :: fileformatX, fileformatY, fileformatZ
    character(len=1),parameter :: NL=char(10) !new line character
    character(len=30) :: filename

    ! Number of columns
    FSX = 1+3+numscalar+3
    FSY = 1+3+numscalar
    FSZ = 1+3+numscalar
    ! All digits or only 6 digits
    if (flag_all_digits) then
       write(fileformatX, '( "(",I4,"(E24.16),A)" )' ) FSX
       write(fileformatY, '( "(",I4,"(E24.16),A)" )' ) FSY
       write(fileformatZ, '( "(",I4,"(E24.16),A)" )' ) FSZ
       FSX = 24*FSX+1
       FSY = 24*FSY+1
       FSZ = 24*FSZ+1
    else
       write(fileformatX, '( "(",I4,"(E14.6),A)" )' ) FSX
       write(fileformatY, '( "(",I4,"(E14.6),A)" )' ) FSY
       write(fileformatZ, '( "(",I4,"(E14.6),A)" )' ) FSZ
       FSX = 14*FSX+1
       FSY = 14*FSY+1
       FSZ = 14*FSZ+1
    endif

    do i = 1, nprobes
      if (rankprobes(i)) then
        write(filename,"('./probes/probe_dx_',I4.4)") i
        open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
             ,access='direct',recl=FSX)
        write(iounit,fileformatX,rec=itime-extra_probes_offset) t, dfdx(i,:), NL
        close(iounit)
      endif
      if (rankprobesY(i)) then
        write(filename,"('./probes/probe_dy_',I4.4)") i
        open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
             ,access='direct',recl=FSY)
        write(iounit,fileformatY,rec=itime-extra_probes_offset) t, dfdy(i,:), NL
        close(iounit)
      endif
      if (rankprobesZ(i)) then
        write(filename,"('./probes/probe_dz_',I4.4)") i
        open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
             ,access='direct',recl=FSZ)
        write(iounit,fileformatZ,rec=itime-extra_probes_offset) t, dfdz(i,:), NL
        close(iounit)
      endif
    enddo

  end subroutine write_extra_probes

  !############################################################################
  subroutine log_probes_position()

    use param, only : xlx, zlz, dx, dz, half
    use variables, only : nxm, yp, ypi, nzm

    integer :: iounit, i
    real(mytype) :: x, y, z, xpp, ypp, zpp

    open(newunit=iounit,file="probes/location.dat")
    write(iounit, *) "#"
    write(iounit, *) "# Location of the probes on the velocity and pressure grids"
    write(iounit, *) "#"
    write(iounit, *) "# Probe Id, xvel, yvel, zvel, xp, yp, zp"
    write(iounit, *) "#"
    do i = 1, nprobes
       x = real(nxprobes(i) - 1, mytype) * dx
       xpp = real(nxprobesP(i) - half, mytype) * xlx / real(nxm, mytype)
       y = yp(nyprobes(i))
       ypp = ypi(nyprobesP(i))
       z = real(nzprobes(i) - 1, mytype) * dz
       zpp = real(nzprobesP(i) - half, mytype) * zlz / real(nzm, mytype)
       write(iounit, *) i, x, y, z, xpp, ypp, zpp
    enddo
    close(iounit)

  end subroutine log_probes_position

  !############################################################################
  !
  ! Free allocated memory at the end of the simulation
  !
  subroutine finalize_probes()

    if (nprobes.le.0) return

    deallocate(xyzprobes)
    deallocate(nxprobes, nyprobes, nzprobes)
    deallocate(rankprobes, rankprobesY, rankprobesZ)
    deallocate(nxprobesP, nyprobesP, nzprobesP)
    deallocate(rankprobesP)
    if (flag_extra_probes) then
      deallocate(dfdx)
      deallocate(dfdy)
      deallocate(dfdz)
    endif

  end subroutine finalize_probes

end module probes
