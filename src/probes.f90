!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: probes.f90
!!!      AUTHOR: Felipe Schuch
!!!    MODIFIED: Cédric Flageul
!!! DESCRIPTION: This module is dedicated to monitoring points.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module probes

  USE decomp_2d, only : ph1, nrank, mytype, xstart, xend

  IMPLICIT NONE

  ! Number of probes
  integer, save :: nprobes
  ! ???
!  integer, save :: ntimes1, ntimes2
  ! Probes location on the velocity grid
  logical, save, allocatable, dimension(:) :: rankprobes
  integer, save, allocatable, dimension(:) :: nxprobes, nyprobes, nzprobes
  ! Probes location on the pressure grid
  logical, save, allocatable, dimension(:) :: rankprobesP
  integer, save, allocatable, dimension(:) :: nxprobesP, nyprobesP, nzprobesP
  ! ???
!  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_probes, write_probes

contains

  !############################################################################
  subroutine init_probes(ep1)

    USE MPI
    USE param, only : dx, dy, dz, nclx, ncly, nclz, xlx, yly, zlz, istret, one, half
    USE variables, only : nx, nxm, ny, yp, nym, ypi, nz, nzm

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1

    logical :: fexists
    double precision :: xprobes, yprobes, zprobes, xyzprobes(3)
    integer :: iounit
    integer :: i,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_probes start'
#endif

! ???
!    ntimes1 = 0
!    ntimes2 = 0

! ???
!    allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
!    allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
!    allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))
!    usum=zero;vsum=zero;wsum=zero
!    uusum=zero;uvsum=zero;uwsum=zero
!    vvsum=zero;vwsum=zero;wwsum=zero

    ! Master rank reads number of probes
    if (nrank.eq.0) then
       inquire(file='probes.prm', exist=fexists)
       if (fexists) then
          open (newunit=iounit,file='probes.prm',status='unknown',form='formatted')
          read (iounit,*) a
          read (iounit,*) nprobes
       else
          nprobes = 0
       endif
    endif
    ! Broadcast
    call MPI_BCAST(nprobes,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)

    ! Exit if no file or no probes
    if (nprobes.le.0) return

    ! Probes on the velocity grid
    allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
    rankprobes(:) = .false.
    ! Probes on the pressure grid
    allocate(nxprobesP(nprobes), nyprobesP(nprobes), nzprobesP(nprobes), rankprobesP(nprobes))
    rankprobesP(:) = .false.

    do i = 1, nprobes
       ! Master rank reads position of probe i
       if (nrank.eq.0) then
          read (iounit,*) xprobes, yprobes, zprobes
          xyzprobes = (/xprobes, yprobes, zprobes/)
       end if
       ! Broadcast
       call MPI_BCAST(xyzprobes,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
       ! Enforce bounds
       xprobes = max(epsilon(xlx), min(xlx*(one-epsilon(xlx)), xyzprobes(1)))
       yprobes = max(epsilon(xlx), min(yly*(one-epsilon(xlx)), xyzprobes(2)))
       zprobes = max(epsilon(xlx), min(zlz*(one-epsilon(xlx)), xyzprobes(3)))

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
       if       (ph1%zst(1) .le. nxprobesP(i) .and. nxprobesP(i) .le. ph1%zen(1)) then
          if    (ph1%zst(2) .le. nyprobesP(i) .and. nyprobesP(i) .le. ph1%zen(2)) then
             if (ph1%zst(3) .le. nzprobesP(i) .and. nzprobesP(i) .le. ph1%zen(3)) then
                rankprobesP(i) = .true.
             endif
          endif
       endif
    enddo

    ! Master rank closes file
    if (nrank.eq.0) close(iounit)

    ! Log the real position of the probes in a file
    if (nrank.eq.0) call log_probes_position()

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_probes ok'
#endif

  end subroutine init_probes

  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,pp3,phi1)

    use param, only : itime, t
    use param, only :  npress
    use var, only : nzmsize
    use variables, only : numscalar

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype), intent(in),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

    integer :: iounit, i, FS, FSP
    character(len=100) :: fileformat, fileformatP
    character(len=1),parameter :: NL=char(10) !new line character
    character(len=30) :: filename

    if (nprobes.le.0) return

    ! Number of columns
    FS = 1+3+numscalar
    FSP = 1+1 ! Pressure grid
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    write(fileformatP, '( "(",I4,"(E14.6),A)" )' ) FSP
    ! Line width
    FS = FS*14+1
    FSP = FSP*14+1

    do i=1, nprobes
       if (rankprobes(i)) then
          write(filename,"('./probes/probe',I4.4)") i
          open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          write(iounit,fileformat,rec=itime) t,&                         !1
               ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
               uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
               uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !numscalar
               NL                                                    !+1
          close(iounit)
       endif
       if (rankprobesP(i)) then
          write(filename,"('./probes/probeP',I4.4)") i
          open(newunit=iounit,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FSP)
          write(iounit,fileformatP,rec=itime) t,&                         !1
               pp3(nxprobesP(i),nyprobesP(i),nzprobesP(i),1),&            !2
               NL                                                    !+1
          close(iounit)
       endif
    enddo

  end subroutine write_probes

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

end module probes
