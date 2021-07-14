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

module cavity

  use decomp_2d
  use variables
  use param

  implicit none

  integer, save :: iocavity

  private ! All functions/subroutines private by default
  public :: init_cavity, boundary_conditions_cavity, &
            postprocess_cavity, momentum_forcing_cavity, &
            finalize_cavity

contains

  !
  ! Initialize velocity and temperature
  !
  subroutine init_cavity (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use MPI

    implicit none

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    ! Local variables
    logical :: read_from_file
    integer :: k, j, i, is, code

    !
    ! Simulation can start from a 3D snapshot
    !
    if (nrank.eq.0) inquire(file="cavity_init_ux", exist=read_from_file)
    call MPI_BCAST(read_from_file,1,MPI_LOGICAL,0,MPI_COMM_WORLD,code)
    if (code.ne.0) call decomp_2d_abort(code, "MPI_BCAST")
    if (read_from_file) then

       if (nrank.eq.0) print *, "Cavity: init from snapshot."
       call decomp_2d_read_one(1,ux1,"cavity_init_ux")
       call decomp_2d_read_one(1,uy1,"cavity_init_uy")
       call decomp_2d_read_one(1,uz1,"cavity_init_uz")
       if (iscalar.eq.1) then
          call decomp_2d_read_one(1,phi1(:,:,:,1),"cavity_init_t")
          if (numscalar.gt.1) then
             phi1(:,:,:,2:numscalar) = zero
          endif
       endif

    else
            
       !
       ! Initial condition: linear temperature profile
       !
       if (iscalar==1) then
          do is = 1, numscalar
             do k = 1, xsize(3)
                do j = 1, xsize(2)
                   do i = 1, xsize(1)
                      !
                      ! 0 <= x <=1 and -0.5 <= T <= 0.5
                      ! x = 0 => T =  0.5
                      ! x = 1 => T = -0.5
                      !
                      phi1(i,j,k,is) = zpfive - (i-1+xstart(1)-1) * dx / xlx
                   enddo
                enddo
             enddo
          enddo
       endif

       !
       ! Initial condition: no velocity
       !
       ux1=zero
       uy1=zero
       uz1=zero

    endif

    return
  end subroutine init_cavity

  !
  ! Apply temperature boundary conditions
  !
  subroutine boundary_conditions_cavity (ux,uy,uz,phi)

    implicit none

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    ! Local variables
    integer :: i, j, k, is

    !
    ! Scalar boundary conditions
    !
    if (iscalar.ne.1) return

    !
    ! Imposed temperature at x = 0
    !
    IF (nclxS1.EQ.2) THEN
       if (xstart(1).eq.1) then
          phi(1,:,:,:) = zpfive
       endif
    ENDIF
    !
    ! Imposed temperature at x = Lx
    !
    IF (nclxSn.EQ.2) THEN
       if (xend(1).eq.nx) then
          phi(xsize(1),:,:,:) = - zpfive 
       endif
    ENDIF

    !
    ! Imposed linear profile at y = 0 and y = Ly
    !
    ! B.C. at y = 0. This is not used when nclyS1 = 1
    !
    IF (nclyS1.EQ.2) THEN
       if (xstart(2).eq.1) then
          do is = 1, numscalar
             do k = 1, xsize(3)
                j = 1
                do i = 1, xsize(1)
                   phi(i,j,k,is) = zpfive - (i-1+xstart(1)-1) * dx / xlx
                enddo
             enddo
          enddo
       endif
    ENDIF
    !
    ! B.C. at y = Ly. This is not used when nclySn = 1
    !
    IF (nclySn.EQ.2) THEN
       if (xend(2).eq.ny) then
          do is = 1, numscalar
             do k = 1, xsize(3)
                j = xsize(2)
                do i = 1, xsize(1)
                   phi(i,j,k,is) = zpfive - (i-1+xstart(1)-1) * dx / xlx
                enddo
             enddo
          enddo
       endif
    ENDIF

    !
    ! Periodicity in Z direction
    !

  end subroutine boundary_conditions_cavity

  !
  ! Extract min / max of various quantities at each time step
  !
  subroutine postprocess_cavity(ux,uy,uz,pp,phi)

    use MPI
    use var, only : nzmsize

    implicit none

    ! Arguments
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype),intent(in),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    ! Local variables
    real(mytype) :: phimax, phimin, umax, umin, vmax, vmin, phim, phi2m
    real(mytype) :: array(6), arrayin(6)
    integer :: code

    ! Compute min & max
    umax = maxval(ux(:,:,:))
    umin =-minval(ux(:,:,:))
    vmax = maxval(uy(:,:,:))
    vmin =-minval(uy(:,:,:))
    phimax = maxval(phi(:,:,:,1))
    phimin =-minval(phi(:,:,:,1))
    arrayin(:) =  (/umin, umax, vmin, vmax, phimin, phimax /)

    ! Compute min & max over the domain
    call MPI_REDUCE(arrayin,array,6,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code.ne.0) call decomp_2d_abort(code, "MPI_REDUCE")

    ! Compute average over the domain
    if (iscalar.eq.1) then

      ! Temperature
      phimax = sum(phi(:,:,:,1))
      call MPI_REDUCE(phimax,phim,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
      if (code.ne.0) call decomp_2d_abort(code, "MPI_REDUCE")

      ! Squared temperature
      phimax = sum(phi(:,:,:,1)*phi(:,:,:,1))
      call MPI_REDUCE(phimax,phi2m,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
      if (code.ne.0) call decomp_2d_abort(code, "MPI_REDUCE")

    else

      ! We should not be here
      phim = zero
      phi2m = zero

    endif

    ! Post-processing + output on master rank
    if (nrank.eq.0) then

      phim = phim / real(nx*ny*nz)
      phi2m = phi2m / real(nx*ny*nz)
      phi2m = phi2m - phim*phim

      ! Print header at first time step
      if (itime.eq.ifirst) then
         open(newunit=iocavity, file='cavity.dat', form='formatted')
         write(iocavity,*) "# umin, umax, vmin, vmax, phimin, phimax, average_T, variance_T"
      endif
      ! Use format to avoid printing all the digits
      write(iocavity, "(8E13.5)") -array(1), array(2), &
                                  -array(3), array(4), &
                                  -array(5), array(6), &
                                  phim, phi2m

    endif

  end subroutine postprocess_cavity

  !
  ! Add the gravity term in the momentum equation
  !
  subroutine momentum_forcing_cavity(duy, phi)

    implicit none

    ! Arguments
    real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3), ntime) :: duy
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if (iscalar.eq.1) duy(:,:,:,1) = duy(:,:,:,1) + phi(:,:,:,1)

  end subroutine momentum_forcing_cavity

  !
  ! Finalize Cavity
  !
  subroutine finalize_cavity()

    implicit none

    if (nrank.eq.0) close(iocavity)

  end subroutine finalize_cavity

end module cavity
