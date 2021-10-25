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

module thomas

  use decomp_2d, only : mytype
  use param, only : one

  implicit none

  private
  public :: xthomas, ythomas, zthomas

  interface xthomas
    module procedure xthomas_0
    module procedure xthomas_12
  end interface xthomas

  interface ythomas
    module procedure ythomas_0
    module procedure ythomas_12
  end interface ythomas

  interface zthomas
    module procedure zthomas_0
    module procedure zthomas_12
  end interface zthomas


contains

  ! Thomas algorithm in X direction (periodicity)
  pure subroutine xthomas_0(tt, rr, ss, ff, fs, fw, alfa, nx, ny, nz)

    implicit none

    real(mytype), intent(inout), dimension(nx,ny,nz) :: tt, rr
    real(mytype), intent(out), dimension(ny,nz) :: ss
    real(mytype), intent(in), dimension(nx):: ff, fs, fw
    real(mytype), intent(in) :: alfa
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k

    call xthomas_12(tt, ff, fs, fw, nx, ny, nz)
    call xthomas_12(rr, ff, fs, fw, nx, ny, nz)
    do concurrent (k=1:nz, j=1:ny)
       ss(j,k) = (    tt(1,j,k)-alfa*tt(nx,j,k)) &
               / (one+rr(1,j,k)-alfa*rr(nx,j,k))
       do concurrent (i=1:nx)
          tt(i,j,k) = tt(i,j,k) - ss(j,k)*rr(i,j,k)
       enddo
    enddo

  end subroutine xthomas_0

  ! Thomas algorithm in X direction
  pure subroutine xthomas_12(tt, ff, fs, fw, nx, ny, nz)

    implicit none

    real(mytype), intent(inout), dimension(nx,ny,nz) :: tt
    real(mytype), intent(in), dimension(nx):: ff, fs, fw
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k

    do concurrent (k=1:nz, j=1:ny)
       do i = 2, nx
          tt(i,j,k) = tt(i,j,k) - tt(i-1,j,k)*fs(i)
       enddo
       tt(nx,j,k) = tt(nx,j,k) * fw(nx)
       do i=nx-1,1,-1
          tt(i,j,k) = (tt(i,j,k)-ff(i)*tt(i+1,j,k)) * fw(i)
       enddo
    enddo

  end subroutine xthomas_12

  ! Thomas algorithm in Y direction (periodicity)
  pure subroutine ythomas_0(tt, rr, ss, ff, fs, fw, alfa, nx, ny, nz)

    implicit none

    real(mytype), intent(inout), dimension(nx,ny,nz) :: tt, rr
    real(mytype), intent(out), dimension(nx,nz) :: ss
    real(mytype), intent(in), dimension(ny):: ff, fs, fw
    real(mytype), intent(in) :: alfa
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k

    call ythomas_12(tt, ff, fs, fw, nx, ny, nz)
    call ythomas_12(rr, ff, fs, fw, nx, ny, nz)
    do concurrent (k=1:nz)
       do concurrent (i=1:nx)
          ss(i,k) = (    tt(i,1,k)-alfa*tt(i,ny,k)) &
                  / (one+rr(i,1,k)-alfa*rr(i,ny,k))
       enddo
       do concurrent (j=1:ny, i=1:nx)
          tt(i,j,k) = tt(i,j,k) - ss(i,k)*rr(i,j,k)
       enddo
    enddo

  end subroutine ythomas_0

  ! Thomas algorithm in Y direction
  pure subroutine ythomas_12(tt, ff, fs, fw, nx, ny, nz)

    implicit none

    real(mytype), intent(inout), dimension(nx,ny,nz) :: tt
    real(mytype), intent(in), dimension(ny):: ff, fs, fw
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k

    do concurrent (k=1:nz)
       do j=2,ny
          do concurrent (i=1:nx)
             tt(i,j,k) = tt(i,j,k) - tt(i,j-1,k)*fs(j)
          enddo
       enddo
       do concurrent (i=1:nx)
          tt(i,ny,k) = tt(i,ny,k) * fw(ny)
       enddo
       do j=ny-1,1,-1
          do concurrent (i=1:nx)
             tt(i,j,k) = (tt(i,j,k)-ff(j)*tt(i,j+1,k)) * fw(j)
          enddo
       enddo
    enddo

  end subroutine ythomas_12

  ! Thomas algorithm in Z direction (periodicity)
  pure subroutine zthomas_0(tt, rr, ss, ff, fs, fw, alfa, nx, ny, nz)

    implicit none

    real(mytype), intent(inout), dimension(nx,ny,nz) :: tt, rr
    real(mytype), intent(out), dimension(nx,ny) :: ss
    real(mytype), intent(in), dimension(nz):: ff, fs, fw
    real(mytype), intent(in) :: alfa
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k

    call zthomas_12(tt, ff, fs, fw, nx, ny, nz)
    call zthomas_12(rr, ff, fs, fw, nx, ny, nz)
    do concurrent (j=1:ny, i=1:nx)
       ss(i,j) = (    tt(i,j,1)-alfa*tt(i,j,nz)) &
               / (one+rr(i,j,1)-alfa*rr(i,j,nz))
    enddo
    do concurrent (k=1:nz, j=1:ny, i=1:nx)
       tt(i,j,k) = tt(i,j,k) - ss(i,j)*rr(i,j,k)
    enddo

  end subroutine zthomas_0

  ! Thomas algorithm in Z direction
  pure subroutine zthomas_12(tt, ff, fs, fw, nx, ny, nz)

    implicit none

    real(mytype), intent(inout), dimension(nx,ny,nz) :: tt
    real(mytype), intent(in), dimension(nz):: ff, fs, fw
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k
       
    do k=2,nz
       do concurrent (j=1:ny, i=1:nx)
          tt(i,j,k) = tt(i,j,k) - tt(i,j,k-1)*fs(k)
       enddo
    enddo 
    do concurrent (j=1:ny, i=1:nx)
       tt(i,j,nz) = tt(i,j,nz) * fw(nz)
    enddo
    do k=nz-1,1,-1
       do concurrent (j=1:ny, i=1:nx)
          tt(i,j,k) = (tt(i,j,k)-ff(k)*tt(i,j,k+1)) * fw(k)
       enddo
    enddo 

  end subroutine zthomas_12

end module thomas
