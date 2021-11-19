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

module user_sim

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_user, boundary_conditions_user, postprocess_user, visu_user

contains

  subroutine init_user (ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d_io
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,is,code
    integer, dimension (:), allocatable :: seed
    integer ::  isize

    if (iscalar==1) then

       phi1(:,:,:,:) = zero
    endif

    if (iin.eq.0) then !empty domain

       if (nrank==0) write(*,*) "Empty initial domain!"

       ux1=zero; uy1=zero; uz1=zero

    endif

    if (iin.eq.1) then !generation of a random noise

       !INIT FOR G AND U=MEAN FLOW + NOISE
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
                uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
                uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
             enddo
          enddo
       enddo

    endif

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_user

  subroutine boundary_conditions_user (ux,uy,uz,phi,ep)

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    IF (nclx1.EQ.2) THEN
    ENDIF
    IF (nclxn.EQ.2) THEN
    ENDIF

    IF (ncly1.EQ.2) THEN
    ENDIF
    IF (nclyn.EQ.2) THEN
    ENDIF

    IF (nclz1.EQ.2) THEN
    ENDIF
    IF (nclzn.EQ.2) THEN
    ENDIF

  end subroutine boundary_conditions_user

  subroutine postprocess_user(ux1,uy1,uz1,phi1,ep1)

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

  end subroutine postprocess_user

  !############################################################################
  !!
  !!  SUBROUTINE: visu_user
  !!      AUTHOR: CF
  !! DESCRIPTION: Performs case-specific visualization
  !!
  !############################################################################
  subroutine visu_user(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

  end subroutine visu_user

end module user_sim
