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

module ludecomp

  use decomp_2d, only : mytype

  implicit none

  private
  public :: ludecomp7, ludecomp9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Interface for septadiag/ nonadiag LU decomposition used in implicit mode
  ! septadiag is for when isecondder is not equal 5 ('classic' second order schemes)
  ! nonadiag is for when isecondder is equal to 5 ('diffusive' second order schemes)
  ! ludecomp7/9_0 is called when ncly=0 (when the matrix is cyclic)
  ! ludecomp7/9_12 is called when ncly=1 or 2
  !
  interface ludecomp7
     module procedure ludecomp7_12
     module procedure ludecomp7_0
  end interface ludecomp7
  
  interface ludecomp9
     module procedure ludecomp9_12
     module procedure ludecomp9_0
  end interface ludecomp9

contains

!*******************************************************************
!Decomposition LU matrix septadiag
  subroutine ludecomp7_12(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ny)
!
!*******************************************************************
    use decomp_2d, only : mytype
    USE param
    implicit none
    
    integer :: i,j,k,ny
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm

    real(mytype), dimension(ny,ny) :: mata, matl, matu, prod

    ggm=zero
    hhm=zero
    ssm=zero
    vvm=zero
    wwm=zero
    zzm=zero

    ggm(1)=aam(1)
    hhm(1)=bbm(1)
    ssm(1)=ccm(1)

    zzm(2)=ddm(2)/ggm(1)
    ggm(2)=aam(2)-zzm(2)*hhm(1)
    hhm(2)=bbm(2)-zzm(2)*ssm(1)
    ssm(2)=ccm(2)-zzm(2)*rrm(1)

    wwm(3)=eem(3)/ggm(1)
    zzm(3)=(ddm(3)-wwm(3)*hhm(1))/ggm(2)
    ggm(3)=aam(3)-wwm(3)*ssm(1)-zzm(3)*hhm(2)
    hhm(3)=bbm(3)-wwm(3)*rrm(1)-zzm(3)*ssm(2)
    ssm(3)=ccm(3)-zzm(3)*rrm(2)

    do j=4,ny
       vvm(j)=qqm(j)/ggm(j-3)
       wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
       zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
       ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
       hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
       ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
    enddo


  end subroutine ludecomp7_12

!*******************************************************************
!Decomposition LU matrix cyclic septadiag
  subroutine ludecomp7_0(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,ny)
!
!*******************************************************************
    use decomp_2d, only : mytype
    USE param
    
    implicit none

    integer :: i,j,k,ny
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm
    real(mytype),dimension(ny), intent(out) :: l1m,l2m,l3m
    real(mytype),dimension(ny), intent(out) :: u1m,u2m,u3m

    real(mytype), dimension(ny,ny) :: mata, matl, matu, prod

    ! a=diag, b=diag+1, c=diag+2, r=diag+3
    !         d=diag-1, e=diag-2, q=diag-3
    ggm=zero;hhm=zero;ssm=zero ! U part (diag, diag+1, diag+2, diag+3=rrm)
    u1m=zero;u2m=zero;u3m=zero ! U cyclic part (last col, last-1, last-2)
    vvm=zero;wwm=zero;zzm=zero ! L part (diag=1, diag-3, diag-2, diag-1)
    l1m=zero;l2m=zero;l3m=zero ! L cyclic part (last row, last-1, last-2)

    ggm(1)=aam(1)
    hhm(1)=bbm(1)
    ssm(1)=ccm(1)
    u1m(1)=ddm(1)
    u2m(1)=eem(1)
    u3m(1)=qqm(1)
    l1m(1)=bbm(ny)/ggm(1)
    l2m(1)=ccm(ny-1)/ggm(1)
    l3m(1)=rrm(ny-2)/ggm(1)

    zzm(2)=ddm(2)/ggm(1)
    ggm(2)=aam(2)-zzm(2)*hhm(1)
    hhm(2)=bbm(2)-zzm(2)*ssm(1)
    ssm(2)=ccm(2)-zzm(2)*rrm(1)
    u1m(2)=eem(2)-zzm(2)*u1m(1)
    u2m(2)=qqm(2)-zzm(2)*u2m(1)
    u3m(2)=      -zzm(2)*u3m(1)
    l1m(2)=(ccm(ny)  -hhm(1)*l1m(1))/ggm(2)
    l2m(2)=(rrm(ny-1)-hhm(1)*l2m(1))/ggm(2)
    l3m(2)=(         -hhm(1)*l3m(1))/ggm(2)

    wwm(3)=eem(3)/ggm(1)
    zzm(3)=(ddm(3)-wwm(3)*hhm(1))/ggm(2)
    ggm(3)=aam(3)-wwm(3)*ssm(1)-zzm(3)*hhm(2)
    hhm(3)=bbm(3)-wwm(3)*rrm(1)-zzm(3)*ssm(2)
    ssm(3)=ccm(3)              -zzm(3)*rrm(2)
    u1m(3)=qqm(3)-wwm(3)*u1m(1)-zzm(3)*u1m(2)
    u2m(3)=      -wwm(3)*u2m(1)-zzm(3)*u2m(2)
    u3m(3)=      -wwm(3)*u3m(1)-zzm(3)*u3m(2)
    l1m(3)=(rrm(ny)-ssm(1)*l1m(1)-hhm(2)*l1m(2))/ggm(3)
    l2m(3)=(       -ssm(1)*l2m(1)-hhm(2)*l2m(2))/ggm(3)
    l3m(3)=(       -ssm(1)*l3m(1)-hhm(2)*l3m(2))/ggm(3)

    do j=4,ny-3
       vvm(j)=qqm(j)/ggm(j-3)
       wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
       zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
       ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
       hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
       ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
       u1m(j)= - vvm(j)*u1m(j-3) - wwm(j)*u1m(j-2) - zzm(j)*u1m(j-1)
       u2m(j)= - vvm(j)*u2m(j-3) - wwm(j)*u2m(j-2) - zzm(j)*u2m(j-1)
       u3m(j)= - vvm(j)*u3m(j-3) - wwm(j)*u3m(j-2) - zzm(j)*u3m(j-1)
       l1m(j)= - (rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1))/ggm(j)
       l2m(j)= - (rrm(j-3)*l2m(j-3)+ssm(j-2)*l2m(j-2)+hhm(j-1)*l2m(j-1))/ggm(j)
       l3m(j)= - (rrm(j-3)*l3m(j-3)+ssm(j-2)*l3m(j-2)+hhm(j-1)*l3m(j-1))/ggm(j)
    enddo

    do j=ny-2,ny
       vvm(j)=qqm(j)/ggm(j-3)
       wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
       zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
       ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
       hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
       ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
    enddo
    j=ny-2
    u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
         - l3m(j-1)*rrm(j-1)
    u2m(j)= - vvm(j)*u2m(j-3)-wwm(j)*u2m(j-2)-zzm(j)*u2m(j-1) &
         - l3m(j-2)*rrm(j-2) - l3m(j-1)*ssm(j-1)
    u3m(j)= - vvm(j)*u3m(j-3)-wwm(j)*u3m(j-2)-zzm(j)*u3m(j-1) &
         - l3m(j-3)*rrm(j-3) - l3m(j-2)*ssm(j-2) - l3m(j-1)*hhm(j-1)
    do k=1,j-1
       u1m(j)=u1m(j) - u1m(k)*l3m(k)
       u2m(j)=u2m(j) - u2m(k)*l3m(k)
       u3m(j)=u3m(j) - u3m(k)*l3m(k)
    enddo
    l1m(j)= - ( rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1) &
         + vvm(ny)*u3m(j-1)+wwm(ny)*u3m(j) ) / ( ggm(j)+u3m(j) )
    l2m(j)= - (rrm(j-3)*l2m(j-3)+ssm(j-2)*l2m(j-2)+hhm(j-1)*l2m(j-1) &
         + vvm(ny-1)*u3m(j-2)+wwm(ny-1)*u3m(j-1)+zzm(ny-1)*u3m(j) ) / (ggm(j)+u3m(j))
    do k=1,j-1
       l1m(j)=l1m(j) - l1m(k)*u3m(k) / (ggm(j)+u3m(j))
       l2m(j)=l2m(j) - l2m(k)*u3m(k) / (ggm(j)+u3m(j))
    enddo
    j=ny-1
    u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
         - l2m(j-2)*rrm(j-2) - l2m(j-1)*ssm(j-1)
    u2m(j)= - vvm(j)*u2m(j-3)-wwm(j)*u2m(j-2)-zzm(j)*u2m(j-1) &
         - l2m(j-3)*rrm(j-3) - l2m(j-2)*ssm(j-2) - l2m(j-1)*hhm(j-1)
    do k=1,j-1
       u1m(j)=u1m(j) - u1m(k)*l2m(k)
       u2m(j)=u2m(j) - u2m(k)*l2m(k)
    enddo
    l1m(j)= - (rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1) &
         + vvm(ny)*u2m(j-2)+wwm(ny)*u2m(j-1)+zzm(ny)*u2m(j) ) / (ggm(j)+u2m(j));
    do k=1,j-1
       l1m(j)=l1m(j) - l1m(k)*u2m(k) / (ggm(j)+u2m(j));
    enddo
    j=ny
    u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
         - l1m(j-3)*rrm(j-3)-l1m(j-2)*ssm(j-2)-l1m(j-1)*hhm(j-1);
    do k=1,j-1
       u1m(j)=u1m(j) - u1m(k)*l1m(k)
    enddo

  end subroutine ludecomp7_0

!*******************************************************************
!Decomposition LU matrix nonadiag
  subroutine ludecomp9_12(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ttm,uum,sssm,zzzm,ny)
!
!*******************************************************************
    use decomp_2d, only : mytype
    use param

    implicit none

    integer :: i,j,k,ny
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm,ttm,uum
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm,sssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm,zzzm

    real(mytype), dimension(ny,ny) :: mata, matl, matu, prod

    ggm=zero;hhm=zero;ssm=zero;sssm=zero
    vvm=zero;wwm=zero;zzm=zero;zzzm=zero

    ggm(1)=aam(1)
    hhm(1)=bbm(1)
    ssm(1)=ccm(1)
    sssm(1)=rrm(1)


    zzzm(2)=ddm(2)/ggm(1)
    ggm(2)=aam(2)-zzzm(2)*hhm(1)
    hhm(2)=bbm(2)-zzzm(2)*ssm(1)
    ssm(2)=ccm(2)-zzzm(2)*sssm(1)
    sssm(2)=rrm(2)-zzzm(2)*ttm(1)

    zzm(3)=eem(3)/ggm(1)
    zzzm(3)=(ddm(3)-zzm(3)*hhm(1))/ggm(2)
    ggm(3)=aam(3)-zzm(3)*ssm(1)-zzzm(3)*hhm(2)
    hhm(3)=bbm(3)-zzm(3)*sssm(1)-zzzm(3)*ssm(2)
    ssm(3)=ccm(3)-zzm(3)*ttm(1)-zzzm(3)*sssm(2)
    sssm(3)=rrm(3)-zzzm(3)*ttm(2)

    wwm(4)=qqm(4)/ggm(1)
    zzm(4)=(eem(4)-wwm(4)*hhm(1))/ggm(2)
    zzzm(4)=(ddm(4)-zzm(4)*hhm(2)-wwm(4)*ssm(1))/ggm(3)

    ggm(4)=aam(4)-wwm(4)*sssm(1)-zzm(4)*ssm(2)-zzzm(4)*hhm(3)
    hhm(4)=bbm(4)-wwm(4)*ttm(1)-zzm(4)*sssm(2)-zzzm(4)*ssm(3)
    ssm(4)=ccm(4)-zzm(4)*ttm(2)-zzzm(4)*sssm(3)
    sssm(4)=rrm(4)-zzzm(4)*ttm(3)

    do j=5,ny
       vvm(j)= uum(j)/ggm(j-4)
       wwm(j)=(qqm(j)-vvm(j)*hhm(j-4))/ggm(j-3)
       zzm(j)=(eem(j)-wwm(j)*hhm(j-3)-vvm(j)*ssm(j-4))/ggm(j-2)
       zzzm(j)=(ddm(j)-zzm(j)*hhm(j-2)-wwm(j)*ssm(j-3)-vvm(j)*sssm(j-4))/ggm(j-1)

       ggm(j)=aam(j)-vvm(j)*ttm(j-4)-wwm(j)*sssm(j-3)-zzm(j)*ssm(j-2)-zzzm(j)*hhm(j-1)
       hhm(j)=bbm(j)-wwm(j)*ttm(j-3)-zzm(j)*sssm(j-2)-zzzm(j)*ssm(j-1)
       ssm(j)=ccm(j)-zzm(j)*ttm(j-2)-zzzm(j)*sssm(j-1)
       sssm(j)=rrm(j)-zzzm(j)*ttm(j-1)

    enddo
  end subroutine ludecomp9_12

!*******************************************************************
!Decomposition LU matrix cyclic nonadiag
  subroutine ludecomp9_0(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,ny)
!
!*******************************************************************
    use decomp_2d, only : mytype
    use param
    USE MPI

    implicit none

    integer :: i,j,k,ny
    integer :: code,ierror
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm
    real(mytype),dimension(ny), intent(out) :: l1m,l2m,l3m
    real(mytype),dimension(ny), intent(out) :: u1m,u2m,u3m

    real(mytype), dimension(ny,ny) :: mata, matl, matu, prod


    ggm=zero;hhm=zero;ssm=zero
    u1m=zero;u2m=zero;u3m=zero
    vvm=zero;wwm=zero;zzm=zero
    l1m=zero;l2m=zero;l3m=zero

    print *,'NOT READY YET! SIMULATION IS STOPPED!'
    call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop

  end subroutine ludecomp9_0

end module ludecomp

!*******************************************************************
!
module matinv

  use decomp_2d, only : mytype

  implicit none

  private
  public :: septinv, nonainv
  !
  ! Interface for septadiag/nonadiag matrix inversion used in implicit mode
  ! septadiag is for when isecondder is not equal 5 ('classic' second order schemes)
  ! nonadiag is for when isecondder is equal to 5 ('diffusive' second order schemes)
  ! septinv_0 is called when ncly=0 (when the matrix is cyclic)
  ! septinv_12 is called when ncly=1 or 2
  !
  interface septinv
     module procedure septinv_12
     module procedure septinv_0
  end interface septinv

  interface nonainv
     module procedure nonainv_12
     module procedure nonainv_0
  end interface nonainv

contains

!********************************************************************
!INVERSE SEPTADIAG MATRIX
  subroutine septinv_12(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,nx,ny,nz)
    !
    !********************************************************************
    use decomp_2d, only : mytype

    implicit none

    integer :: i,j,k,nx,ny,nz
    real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
    real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
    real(mytype),dimension(ny),intent(in)    :: ggm,hhm,ssm,rrm
    real(mytype),dimension(ny),intent(in)    :: vvm,wwm,zzm
    
    do k=1,nz
       do i=1,nx

          !going down
          xsol(i,1,k)=bbb(i,1,k)
          xsol(i,2,k)=bbb(i,2,k)-zzm(2)*xsol(i,1,k)
          xsol(i,3,k)=bbb(i,3,k)-zzm(3)*xsol(i,2,k)-wwm(3)*xsol(i,1,k)

          do j=4,ny
             xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                  -zzm(j)*xsol(i,j-1,k)
          enddo
          !

          !going up
          xsol(i,ny,k)=xsol(i,ny,k)/ggm(ny)
          xsol(i,ny-1,k)=(xsol(i,ny-1,k)-hhm(ny-1)*xsol(i,ny,k))/ggm(ny-1)
          xsol(i,ny-2,k)=(xsol(i,ny-2,k)-hhm(ny-2)*xsol(i,ny-1,k)- &
               ssm(ny-2)*xsol(i,ny,k))/ggm(ny-2)

          do j=ny-3,1,-1
             xsol(i,j,k)=(xsol(i,j,k)-hhm(j)*xsol(i,j+1,k)-ssm(j)*xsol(i,j+2,k) &
                  -rrm(j)*xsol(i,j+3,k))/ggm(j)
          enddo

       enddo
    enddo

  end subroutine septinv_12

!********************************************************************
!INVERSE SEPTADIAG CYCLIC
  subroutine septinv_0(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,nx,ny,nz)
!    
!********************************************************************
    use decomp_2d, only : mytype

    implicit none

    integer :: i,j,k,kk,nx,ny,nz
    real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
    real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
    real(mytype),dimension(ny), intent(in)    :: ggm,hhm,ssm,rrm
    real(mytype),dimension(ny), intent(in)    :: vvm,wwm,zzm
    real(mytype),dimension(ny), intent(in) :: l1m,l2m,l3m
    real(mytype),dimension(ny), intent(in) :: u1m,u2m,u3m

    do k=1,nz
       do i=1,nx

          !going down with triangle matrix inf
          xsol(i,1,k)=bbb(i,1,k)
          xsol(i,2,k)=bbb(i,2,k)-zzm(2)*xsol(i,1,k)
          xsol(i,3,k)=bbb(i,3,k)-zzm(3)*xsol(i,2,k)-wwm(3)*xsol(i,1,k)
          do j=4,ny-3
             xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                  -zzm(j)*xsol(i,j-1,k)
          enddo
          j=ny-2
          xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
               -zzm(j)*xsol(i,j-1,k)
          do kk=1,j-1
             xsol(i,j,k)=xsol(i,j,k)-l3m(kk)*xsol(i,kk,k)
          enddo
          j=ny-1
          xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
               -zzm(j)*xsol(i,j-1,k)
          do kk=1,j-1
             xsol(i,j,k)=xsol(i,j,k)-l2m(kk)*xsol(i,kk,k)
          enddo
          j=ny
          xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
               -zzm(j)*xsol(i,j-1,k)
          do kk=1,j-1
             xsol(i,j,k)=xsol(i,j,k)-l1m(kk)*xsol(i,kk,k)
          enddo
          !

          !going up with triangle matrix up
          xsol(i,ny,k)=xsol(i,ny,k)/(ggm(ny)+u1m(ny))
          xsol(i,ny-1,k)=(xsol(i,ny-1,k)-(hhm(ny-1)+u1m(ny-1))*xsol(i,ny,k))/(ggm(ny-1)+u2m(ny-1))
          xsol(i,ny-2,k)=(xsol(i,ny-2,k)-(hhm(ny-2)+u2m(ny-2))*xsol(i,ny-1,k)- &
               (ssm(ny-2)+u1m(ny-2))*xsol(i,ny,k))/(ggm(ny-2)+u3m(ny-2))

          do j=ny-3,1,-1
             xsol(i,j,k)=(xsol(i,j,k)-hhm(j)*xsol(i,j+1,k)-ssm(j)*xsol(i,j+2,k) &
                  -rrm(j)*xsol(i,j+3,k) &
                  -u3m(j)*xsol(i,ny-2,k) &
                  -u2m(j)*xsol(i,ny-1,k) &
                  -u1m(j)*xsol(i,ny,k) )/ggm(j)
          enddo
       enddo
    enddo

  end subroutine septinv_0
  
!********************************************************************
!INVERSE NONADIAG
  subroutine nonainv_12(xSol,bbb,ggm,hhm,ssm,sssm,ttm,zzzm,zzm,wwm,vvm,nx,ny,nz)
!
!********************************************************************
    use decomp_2d, only : mytype

    implicit none

    integer :: i,j,k,nx,ny,nz
    real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
    real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
    real(mytype),dimension(ny),intent(in)    :: ggm,hhm,ssm, sssm, ttm
    real(mytype),dimension(ny),intent(in)    :: vvm,wwm,zzm,zzzm

    do k=1,nz
       do i=1,nx

          !going down
          xSol(i,1,k)=bbb(i,1,k);
          xSol(i,2,k)=bbb(i,2,k)-zzzm(2)*xSol(i,1,k);
          xSol(i,3,k)=bbb(i,3,k)-zzm(3)*xSol(i,1,k)-zzzm(3)*xSol(i,2,k);
          xSol(i,4,k)=bbb(i,4,k)-wwm(4)*xSol(i,1,k)-zzm(4)*xSol(i,2,k)-zzzm(4)*xSol(i,3,k);

          do j=5,ny
             xSol(i,j,k)=bbb(i,j,k)-vvm(j)*xSol(i,j-4,k)-wwm(j)*xSol(i,j-3,k) &
                  -zzm(j)*xSol(i,j-2,k)-zzzm(j)*xSol(i,j-1,k);
          enddo
          !

          !going up
          xSol(i,ny,k)=xSol(i,ny,k)/ggm(ny);
          xSol(i,ny-1,k)=(xSol(i,ny-1,k)-hhm(ny-1)*xSol(i,ny,k))/ggm(ny-1);
          xSol(i,ny-2,k)=(xSol(i,ny-2,k)-hhm(ny-2)*xSol(i,ny-1,k)-ssm(ny-2)*xSol(i,ny,k))/ggm(ny-2);
          xSol(i,ny-3,k)=(xSol(i,ny-3,k)-hhm(ny-3)*xSol(i,ny-2,k)-ssm(ny-3)*xSol(i,ny-1,k) &
               -sssm(ny-3)*xSol(i,ny,k))/ggm(ny-3);

          do j=ny-4,1,-1
             xSol(i,j,k)=(xSol(i,j,k)-hhm(j)*xSol(i,j+1,k)-ssm(j)*xSol(i,j+2,k) &
                  -sssm(j)*xSol(i,j+3,k)-ttm(j)*xSol(i,j+4,k))/ggm(j);
          enddo

       enddo
    enddo

  end subroutine nonainv_12

!********************************************************************
!INVERSE NONADIAG CYCLIC
  subroutine nonainv_0(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,nx,ny,nz)
    !
    !********************************************************************
    use decomp_2d, only : mytype
    USE MPI
    
    implicit none

    integer :: i,j,k,kk,nx,ny,nz
    integer :: code,ierror
    real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
    real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
    real(mytype),dimension(ny), intent(in)    :: ggm,hhm,ssm,rrm
    real(mytype),dimension(ny), intent(in)    :: vvm,wwm,zzm
    real(mytype),dimension(ny), intent(in) :: l1m,l2m,l3m
    real(mytype),dimension(ny), intent(in) :: u1m,u2m,u3m

    print *,'NOT READY YET! SIMULATION IS STOPPED!'
    call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop

  end subroutine nonainv_0
  
end module matinv

!********************************************************************
!
subroutine  inttimp (var1,dvar1,forcing1)
  !
  !********************************************************************
  USE param
  USE variables
  USE var, ONLY: ta1, ta2, tb2, tc2, td2
  USE decomp_2d
  use derivY
  use matinv

  implicit none

  integer :: i,j,k,code

  !! IN
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: forcing1

  !! IN/OUT
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1
  real(mytype),dimension(ysize(1),ysize(3)) :: bc1


!!!!!!!!!!!!!!!
!!! CN2+AB3 !!!
!!!!!!!!!!!!!!!
  if (itimescheme==7) then
     if ((itime.eq.1).and.(irestart.eq.0)) then

        if (nrank==0) print *,'start AB3 with Euler',itime
        !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES

        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 !uhat
                 ta1(i,j,k) = dt*dvar1(i,j,k,1)-forcing1(i,j,k)

                 !save (N+Lxz)(n-1)
                 dvar1(i,j,k,2)=dvar1(i,j,k,1)
              enddo
           enddo
        enddo

     else if ((itime.eq.2).and.(irestart.eq.0)) then

        if (nrank==0) print *,'then with AB2',itime
        !THEN AB2 + IMPLICIT CN SCHEMES
        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 !uhat
                 ta1(i,j,k)= 1.5*dt*dvar1(i,j,k,1)-0.5*dt*dvar1(i,j,k,2)-forcing1(i,j,k)

                 !save (N+Lxz)(n-1)&(n-2)
                 dvar1(i,j,k,3)=dvar1(i,j,k,2)

                 dvar1(i,j,k,2)=dvar1(i,j,k,1)

              enddo
           enddo
        enddo

     else !FINALLY EXPLICIT AB3 + IMPLICIT CN SCHEMES

        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 !uhat
                 ta1(i,j,k)= adt(itr)*dvar1(i,j,k,1)+bdt(itr)*dvar1(i,j,k,2) &
                      +cdt(itr)*dvar1(i,j,k,3)-forcing1(i,j,k)

                 !save (N+Lxz)(n-1)
                 dvar1(i,j,k,3)=dvar1(i,j,k,2)

                 dvar1(i,j,k,2)=dvar1(i,j,k,1)

              enddo
           enddo
        enddo

     endif
  endif

  !Y-PENCIL FOR MATRIX INVERSION
  call transpose_x_to_y(var1,tb2)

  call transpose_x_to_y(ta1,ta2)

  !BC FOR THE VELOCITY (IF NOTHING IS DONE THEN ZERO IS IMPOSED
  if (itype.eq.itype_tbl) then
     bc1(:,:)=tb2(:,ny-1,:)
     !in order to mimick a Neumann BC at the top of the domain for the TBL
  endif
 
  !ta2: A.uhat
  !td2:(A+xcstB).un
  !if isecondder=5, we need nona inversion
  !id isecondder is not 5, we need septa inversion

  !NEED TO BE DONE: CASE 1-2 and 2-1!!!!
  
  if (isecondder.ne.5) then
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        call multmatrix7(td2,ta2,tb2,0)
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        call multmatrix7(td2,ta2,tb2,1)
!        call multmatrix7(te2,tb2,uy2,0) !NEED TO BE DONE: ADJUST PARITY WHEN TB2=UY, NPAIRE in parameters 
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call multmatrix7(td2,ta2,tb2,0)
     endif
  endif
  
  if (isecondder.eq.5) then
     !TO BE DONE: Different types of BC
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        !NOT READY
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        !NOT READY
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call multmatrix9(td2,ta2,tb2)
     endif
  endif

  !Full second member BBB=A uhat+(A+xcst.B)u^n
  ta2(:,:,:)=td2(:,:,:)+ta2(:,:,:)

  !BC for the TBL case
  if (itype.eq.itype_tbl) then
     ta2(:,ny,:)=bc1(:,:)
  endif
 
  !Inversion of the linear system Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
  !if secondder=5, we need nona inversion
  !if isecondder is not 5, we need septa inversion

  !NEED TO BE DONE: CASE 1-2 and 2-1!!!!
  if (isecondder.ne.5) then
     tb2=0.;
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        call septinv(tb2,ta2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        call septinv(tb2,ta2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call septinv(tb2,ta2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
     endif
  endif
  if (isecondder.eq.5) then
     tb2=0.;
     !TO BE DONE: Different types of BC
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        !NOT READY
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        !NOT READY
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call nonainv(tb2,ta2,ggm,hhm,ssm,sssm,ttm,zzzm,zzm,wwm,vvm,ysize(1),ysize(2),ysize(3))
     endif
  endif
  
  !  if (ncly1.eq.0) then
!     call septinv(tb2,ta2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
!  elseif (ncly1.eq.1) then
!     call septinv(tb2,ta2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
!  elseif (ncly1.eq.2) then
!     if (isecondder.ne.5) then
!        call septinv(tb2,ta2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
!     else
!        call nonainv(tb2,ta2,ggm,hhm,ssm,sssm,ttm,zzzm,zzm,wwm,vvm,ysize(1),ysize(2),ysize(3))
!     endif
!  endif

  call transpose_y_to_x(tb2,var1)

  if ( (irestart.eq.1).or.(itime.gt.1) ) then
     var1(:,:,:)=var1(:,:,:)+forcing1(:,:,:)
  endif


  return
end subroutine inttimp

!********************************************************************
!PREMULT FOR INTTIMP WITH SEPTADIAG
subroutine multmatrix7(td2,ta2,ux2,npaire)
! 
!********************************************************************
USE param
USE variables
USE derivY
USE decomp_2d

implicit none

integer :: npaire
integer :: i,j,k,code
!real(mytype) :: xcst ! modification module param
real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2

!xcst= xnu*gdt(itr)*0.5

!A.uhat
   if (istret.ne.0) then
      do j=1,ysize(2)
         ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
      enddo
   endif

   if ((ncly1.eq.0).and.(nclyn.eq.0)) then

      td2(:,1,:) = alsajy*ta2(:,2,:) + ta2(:,1,:) + alsajy*ta2(:,ysize(2),:)
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      td2(:,ysize(2),:) = alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) + alsajy*ta2(:,1,:)
      ta2=td2

   elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then

      if (npaire.eq.0) then
         td2(:,1,:) = ta2(:,1,:)
      else
         td2(:,1,:) = 2.*alsajy*ta2(:,2,:) + ta2(:,1,:)
      endif

      do j=2,ysize(2)-1
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo

      if (npaire.eq.0) then
         td2(:,ysize(2),:) = ta2(:,ysize(2),:)
      else
         td2(:,ysize(2),:) = 2.*alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:)
      endif
      ta2=td2

   elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then

      td2(:,1,:) = 0.
      td2(:,2,:) = alsa2y*ta2(:,1,:) + ta2(:,2,:) + alsa2y*ta2(:,3,:)
      td2(:,3,:) = alsa3y*ta2(:,2,:) + ta2(:,3,:) + alsa3y*ta2(:,4,:)
      do j=4,ysize(2)-3
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      td2(:,ysize(2)-2,:) = alsaty*ta2(:,ysize(2)-3,:) + ta2(:,ysize(2)-2,:) + alsaty*ta2(:,ysize(2)-1,:)
      td2(:,ysize(2)-1,:) = alsamy*ta2(:,ysize(2)-2,:) + ta2(:,ysize(2)-1,:) + alsamy*ta2(:,ysize(2),:)
      td2(:,ysize(2),:) = 0.
      ta2=td2

   endif


!(A+nu*dt.B).un
   if ((ncly1.eq.0).and.(nclyn.eq.0)) then

      call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

   elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then

      if (npaire==0) then
         call deryy(td2,ux2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
      else
         call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
      endif

   elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then

      call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

   endif

   if (istret.ne.0) then
      do j=1,ysize(2)
         ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
      enddo
   endif

   if ((ncly1.eq.0).and.(nclyn.eq.0)) then

      td2(:,1,:) = alsajy*ux2(:,2,:) + ux2(:,1,:) + alsajy*ux2(:,ysize(2),:) &
                    + xcst*td2(:,1,:)
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + xcst*td2(:,j,:)
      enddo
      td2(:,ysize(2),:) = alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajy*ux2(:,1,:) &
                           + xcst*td2(:,ysize(2),:)

   elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then

      if (npaire.eq.0) then
         td2(:,1,:) = ux2(:,1,:) + xcst*td2(:,1,:)
      else
         td2(:,1,:) = 2.*alsajy*ux2(:,2,:) + ux2(:,1,:) &
                    + xcst*td2(:,1,:)
      endif

      do j=2,ysize(2)-1
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + xcst*td2(:,j,:)
      enddo

      if (npaire.eq.0) then
         td2(:,ysize(2),:) = ux2(:,ysize(2),:) + xcst*td2(:,ysize(2),:)
      else
         td2(:,ysize(2),:) = 2.*alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) &
                           + xcst*td2(:,ysize(2),:)
      endif

   elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then

      td2(:,1,:) = 0.
      td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
                 + xcst*td2(:,2,:)
      td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
                 + xcst*td2(:,3,:)
      do j=4,ysize(2)-3
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + xcst*td2(:,j,:)
      enddo
      td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
                          + xcst*td2(:,ysize(2)-2,:)
      td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
                          + xcst*td2(:,ysize(2)-1,:)
      td2(:,ysize(2),:) = 0.

   endif

end subroutine multmatrix7

!********************************************************************

subroutine multmatrix9(td2,ta2,ux2)
  !
  !********************************************************************
  USE param
  USE variables
  USE derivY
  USE decomp_2d

  implicit none

  integer :: i,j,k,code
  !real(mytype) :: xcst ! modification module param
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2


  !A.uhat

  if (istret.ne.0) then
     do j=1,ysize(2)
        ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
     enddo
  endif

  if (ncly1.eq.0) then

     td2(:,1,:) = alsajy*ta2(:,2,:) + ta2(:,1,:) + alsajy*ta2(:,ysize(2),:)
     do j=2,ysize(2)-1
        td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
     enddo
     td2(:,ysize(2),:) = alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) + alsajy*ta2(:,1,:)
     ta2=td2

  elseif (ncly1.eq.1) then


  elseif (ncly1.eq.2) then

     td2(:,1,:) = 0.
     td2(:,2,:) = alsa2y*ta2(:,1,:) + ta2(:,2,:) + alsa2y*ta2(:,3,:)
     td2(:,3,:) = alsa3y*ta2(:,2,:) + ta2(:,3,:) + alsa3y*ta2(:,4,:)
     td2(:,4,:) = alsa4y*ta2(:,3,:) + ta2(:,4,:) + alsa4y*ta2(:,5,:)
     do j=5,ysize(2)-4
        td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
     enddo
     td2(:,ysize(2)-3,:) = alsatty*ta2(:,ysize(2)-4,:) + ta2(:,ysize(2)-3,:) + alsatty*ta2(:,ysize(2)-2,:)
     td2(:,ysize(2)-2,:) = alsaty*ta2(:,ysize(2)-3,:) + ta2(:,ysize(2)-2,:) + alsaty*ta2(:,ysize(2)-1,:)
     td2(:,ysize(2)-1,:) = alsamy*ta2(:,ysize(2)-2,:) + ta2(:,ysize(2)-1,:) + alsamy*ta2(:,ysize(2),:)
     td2(:,ysize(2),:) = 0.
     ta2=td2

  endif

  !(A+nu*dt.B).un

  if (ncly1.eq.0) then

     call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

  elseif (ncly1.eq.1) then


  elseif (ncly1.eq.2) then

     call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

  endif

  if (istret.ne.0) then
     do j=1,ysize(2)
        ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
     enddo
  endif

  if (ncly1.eq.0) then

     td2(:,1,:) = alsajy*ux2(:,2,:) + ux2(:,1,:) + alsajy*ux2(:,ysize(2),:) &
          + xcst*td2(:,1,:)
     do j=2,ysize(2)-1
        td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
             + xcst*td2(:,j,:)
     enddo
     td2(:,ysize(2),:) = alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajy*ux2(:,1,:) &
          + xcst*td2(:,ysize(2),:)

  elseif (ncly1.eq.1) then

  elseif (ncly1.eq.2) then

     td2(:,1,:) = 0.
     td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
          + xcst*td2(:,2,:)
     td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
          + xcst*td2(:,3,:)
     td2(:,4,:) = alsa4y*ux2(:,3,:) + ux2(:,4,:) + alsa4y*ux2(:,5,:) &
          + xcst*td2(:,4,:)
     do j=5,ysize(2)-4
        td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
             + xcst*td2(:,j,:)
     enddo
     td2(:,ysize(2)-3,:) = alsatty*ux2(:,ysize(2)-4,:) + ux2(:,ysize(2)-3,:) + alsatty*ux2(:,ysize(2)-2,:) &
          + xcst*td2(:,ysize(2)-3,:)
     td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
          + xcst*td2(:,ysize(2)-2,:)
     td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
          + xcst*td2(:,ysize(2)-1,:)
     td2(:,ysize(2),:) = 0.

  endif

end subroutine multmatrix9

!********************************************************************
!
subroutine implicit_schemes()
  !
  !********************************************************************

  USE param
  USE derivY
  USE variables
  USE var
  USE param
  use ludecomp

  implicit none

  integer  :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MATRIX M2 TIME IMPLICIT !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!FOR VELOCITY  !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

xcst = xnu*dt*half

if (isecondder.ne.5) then
!!!!!!!!!!!!!!!!!!!!!!
     !7-DIAG !!
!!!!!!!!!!!!!!!!!!!!!!

!!! NCL = 2, dirichlet BC
     !
     !DIAG
     aam(1     )=as1y
     aam(ny    )=asny
     aam(2     )=-two*as2y
     aam(ny-1  )=-two*asmy
     aam(3     )=-two*(as3y+bs3y)
     aam(ny-2  )=-two*(asty+bsty)
     aam(4:ny-3)=-two*(asjy+bsjy+csjy)
     if (istret==0) then
        aam = 1.-xcst*aam
     else
        aam = one/pp2y -xcst*aam
     endif
     !BC on aam
     aam(1 )=one
     aam(ny)=one
     !
     !DIAG SUP 1
     bbm(1     )=bs1y
     bbm(ny    )=bsny
     bbm(2     )=as2y
     bbm(ny-1  )=asmy
     bbm(3     )=as3y
     bbm(ny-2  )=asty
     bbm(4:ny-3)=asjy
     bbm = -xcst*bbm
     if (istret==0) then
        bbm(2     )=bbm(2     )+alsa2y
        bbm(ny-1  )=bbm(ny-1  )+alsamy
        bbm(3     )=bbm(3     )+alsa3y
        bbm(ny-2  )=bbm(ny-2  )+alsaty
        bbm(4:ny-3)=bbm(4:ny-3)+alsajy
     else
        bbm(2     )=bbm(2     )+alsa2y/pp2y(3)
        bbm(ny-1  )=bbm(ny-1  )+alsamy/pp2y(ny)
        bbm(3     )=bbm(3     )+alsa3y/pp2y(4)
        bbm(ny-2  )=bbm(ny-2  )+alsaty/pp2y(ny-1)
        bbm(4:ny-3)=bbm(4:ny-3)+alsajy/pp2y(5:ny-2)
     endif
     !BC on bbm
     bbm(1 )=zero
     bbm(ny)=zero
     !
     !DIAG SUP 2
     ccm(1     )=cs1y
     ccm(ny    )=csny
     ccm(2     )=zero
     ccm(ny-1  )=zero
     ccm(3     )=bs3y
     ccm(ny-2  )=bsty
     ccm(4:ny-3)=bsjy
     ccm = -xcst*ccm
     !BC on ccm
     ccm(1 )=zero
     ccm(ny)=zero
     !
     !DIAG SUP 3
     rrm(1     )=ds1y
     rrm(ny    )=dsny
     rrm(2     )=zero
     rrm(ny-1  )=zero
     rrm(3     )=zero
     rrm(ny-2  )=zero
     rrm(4:ny-3)=csjy
     rrm = -xcst*rrm
     !BC on rrm
     rrm(1 )=zero
     rrm(ny)=zero
     !
     !DIAG INF 1
     if (istret==0) then
        ddm=bbm
     else
        ddm(1     )=bs1y
        ddm(ny    )=bsny
        ddm(2     )=as2y
        ddm(ny-1  )=asmy
        ddm(3     )=as3y
        ddm(ny-2  )=asty
        ddm(4:ny-3)=asjy
        ddm = -xcst*ddm
        ddm(2     )=ddm(2     )+alsa2y/pp2y(1)
        ddm(ny-1  )=ddm(ny-1  )+alsamy/pp2y(ny-2)
        ddm(3     )=ddm(3     )+alsa3y/pp2y(2)
        ddm(ny-2  )=ddm(ny-2  )+alsaty/pp2y(ny-3)
        ddm(4:ny-3)=ddm(4:ny-3)+alsajy/pp2y(3:ny-4)
        !BC on ddm
        ddm(1 )=zero
        ddm(ny)=zero
     endif
     !
     !DIAG INF 2
     eem=ccm
     !
     !DIAG INF 3
     qqm=rrm

!!! NCL = 1, npaire=0, neumann BC, odd function
     !
     ! DIAG
     aam10(1     )=zero
     aam10(ny    )=zero
     aam10(2     )=-two*asjy-three*bsjy-two*csjy
     aam10(ny-1  )=aam10(2)
     aam10(3:ny-2)=-two*(asjy+bsjy+csjy)
     if (istret==0) then
        aam10 = one - xcst*aam10
     else
        aam10 = one/pp2y - xcst*aam10
     endif
     !
     !DIAG SUP 1
     bbm10(1     )=zero
     bbm10(ny    )=zero
     bbm10(2     )=asjy-csjy
     bbm10(ny-1  )=asjy
     bbm10(3     )=asjy
     bbm10(ny-2  )=asjy-csjy
     bbm10(4:ny-3)=asjy
     if (istret==0) then
        bbm10 = alsajy - xcst*bbm10
     else
        bbm10(2:ny-1) = alsajy/pp2y(3:ny) - xcst*bbm10(2:ny-1)
     endif
     !BC on bbm10
     bbm10(1 )=zero
     bbm10(ny)=zero
     !
     !DIAG SUP 2
     ccm10(1     )=zero
     ccm10(ny    )=zero
     ccm10(2     )=bsjy
     ccm10(ny-1  )=zero
     ccm10(3     )=bsjy
     ccm10(ny-2  )=bsjy
     ccm10(4:ny-3)=bsjy
     ccm10 = -xcst*ccm10
     !
     !DIAG SUP 3
     rrm10(1     )=zero
     rrm10(ny    )=zero
     rrm10(2     )=csjy
     rrm10(ny-1  )=zero
     rrm10(3     )=csjy
     rrm10(ny-2  )=zero
     rrm10(4:ny-3)=csjy
     rrm10 = -xcst*rrm10
     !
     !DIAG INF 1
     ddm10(1     )=zero
     ddm10(ny    )=zero
     ddm10(2     )=asjy
     ddm10(ny-1  )=asjy-csjy
     ddm10(3     )=asjy-csjy
     ddm10(ny-2  )=asjy
     ddm10(4:ny-3)=asjy
     if (istret==0) then
        ddm10 = alsajy - xcst*ddm10
     else
        ddm10(2:ny-1) = alsajy/pp2y(1:ny-2) - xcst*ddm10(2:ny-1)
     endif
     !BC on ddm10
     ddm10(1 )=zero
     ddm10(ny)=zero
     !
     !DIAG INF 2
     eem10(1     )=zero
     eem10(ny    )=zero
     eem10(2     )=zero
     eem10(ny-1  )=bsjy
     eem10(3     )=bsjy
     eem10(ny-2  )=bsjy
     eem10(4:ny-3)=bsjy
     eem10 = -xcst*eem10
     !
     !DIAG INF 3
     qqm10(1     )=zero
     qqm10(ny    )=zero
     qqm10(2     )=zero
     qqm10(ny-1  )=csjy
     qqm10(3     )=zero
     qqm10(ny-2  )=csjy
     qqm10(4:ny-3)=csjy
     qqm10 = -xcst*qqm10

!!! NCL = 1, npaire=1, neumann BC, even function
     !
     ! DIAG
     aam11(1     )=-two*(asjy+bsjy+csjy)
     aam11(ny    )=aam11(1)
     aam11(2     )=-two*asjy-bsjy-two*csjy
     aam11(ny-1  )=aam11(2)
     aam11(3:ny-2)=-two*(asjy+bsjy+csjy)
     if (istret==0) then
        aam11 = one - xcst*aam11
     else
        aam11 = one/pp2y - xcst*aam11
     endif
     !
     !DIAG SUP 1
     bbm11(1     )=two*asjy
     bbm11(ny    )=zero
     bbm11(2     )=asjy+csjy
     bbm11(ny-1  )=asjy
     bbm11(3     )=asjy
     bbm11(ny-2  )=asjy+csjy
     bbm11(4:ny-3)=asjy
     if (istret==0) then
        bbm11 = alsajy - xcst*bbm11
     else
        bbm11(1:ny-1) = alsajy/pp2y(2:ny) - xcst*bbm11(1:ny-1)
     endif
     !BC on bbm11
     if (istret==0) then
        bbm11(1 )=bbm11(1)+alsajy
     else
        bbm11(1 )=bbm11(1)+alsajy/pp2y(2)
     endif
     bbm11(ny)=zero
     !
     !DIAG SUP 2
     ccm11(1     )=two*bsjy
     ccm11(ny    )=zero
     ccm11(2     )=bsjy
     ccm11(ny-1  )=zero
     ccm11(3     )=bsjy
     ccm11(ny-2  )=bsjy
     ccm11(4:ny-3)=bsjy
     ccm11 = -xcst*ccm11
     !
     !DIAG SUP 3
     rrm11(1     )=two*csjy
     rrm11(ny    )=zero
     rrm11(2     )=csjy
     rrm11(ny-1  )=zero
     rrm11(3     )=csjy
     rrm11(ny-2  )=zero
     rrm11(4:ny-3)=csjy
     rrm11 = -xcst*rrm11
     !
     !DIAG INF 1
     ddm11(1     )=zero
     ddm11(ny    )=two*asjy
     ddm11(2     )=asjy
     ddm11(ny-1  )=asjy+csjy
     ddm11(3     )=asjy+csjy
     ddm11(ny-2  )=asjy
     ddm11(4:ny-3)=asjy
     if (istret==0) then
        ddm11 = alsajy - xcst*ddm11
     else
        ddm11(2:ny) = alsajy/pp2y(1:ny-1) - xcst*ddm11(2:ny)
     endif
     !BC on ddm11
     ddm11(1 )=zero
     if (istret==0) then
        ddm11(ny)=ddm11(ny)+alsajy!a1
     else
        ddm11(ny)=ddm11(ny)+alsajy/pp2y(ny-1)!a1
     endif
     !
     !DIAG INF 2
     eem11(1     )=zero
     eem11(ny    )=two*bsjy
     eem11(2     )=zero
     eem11(ny-1  )=bsjy
     eem11(3     )=bsjy
     eem11(ny-2  )=bsjy
     eem11(4:ny-3)=bsjy
     eem11 = -xcst*eem11
     !
     !DIAG INF 3
     qqm11(1     )=zero
     qqm11(ny    )=two*csjy
     qqm11(2     )=zero
     qqm11(ny-1  )=csjy
     qqm11(3     )=zero
     qqm11(ny-2  )=csjy
     qqm11(4:ny-3)=csjy
     qqm11 = -xcst*qqm11

!!! NXL = 0
     !DIAG
     if (istret==0) then
        aam0 = one-xcst*(-two*(asjy+bsjy+csjy))
     else
        aam0 = one/pp2y-xcst*(-two*(asjy+bsjy+csjy))
     endif
     !
     !DIAG SUP 1
     if (istret==0) then
        bbm0 = alsajy-xcst*asjy
     else
        bbm0(1:ny-1) = alsajy/pp2y(2:ny) -xcst*asjy
        bbm0(ny) = alsajy/pp2y(1) -xcst*asjy
     endif
     !
     !DIAG SUP 2
     ccm0 = -xcst*bsjy
     !
     !DIAG SUP 3
     rrm0 = -xcst*csjy
     !
     !DIAG INF 1
     if (istret==0) then
        ddm0=bbm0
     else
        ddm0(1) = alsajy/pp2y(ny) -xcst*asjy
        ddm0(2:ny) = alsajy/pp2y(1:ny-1) -xcst*asjy
     endif
     !
     !DIAG INF 2
     eem0=ccm0
     !
     !DIAG INF 3
     qqm0=rrm0

  else
!!!!!!!!!!!!!!!!!!!!!!
     !9-DIAG !!
!!!!!!!!!!!!!!!!!!!!!!

!!! NCL = 2, dirichlet BC
     !
     !DIAG
     aam(1     )=as1y
     aam(ny    )=asny
     aam(2     )=-two*as2y
     aam(ny-1  )=-two*asmy
     aam(3     )=-two*(as3y+bs3y)
     aam(4     )=-two*(as4y+bs4y+cs4y)
     aam(ny-2  )=-two*(asty+bsty)
     aam(ny-3  )=-two*(astty+bstty+cstty)
     aam(5:ny-4)=-two*(asjy+bsjy+csjy+dsjy)
     if (istret==0) then
        aam = one-xcst*aam
     else
        aam = one/pp2y -xcst*aam
     endif
     !BC on aam
     aam(1 )=one
     aam(ny)=one
     !
     !DIAG SUP 1
     bbm(1     )=bs1y
     bbm(ny    )=bsny
     bbm(2     )=as2y
     bbm(ny-1  )=asmy
     bbm(3     )=as3y
     bbm(ny-2  )=asty
     bbm(4     )=as4y
     bbm(ny-3  )=astty
     bbm(5:ny-4)=asjy
     bbm = -xcst*bbm
     if (istret==0) then
        bbm(2     )=bbm(2     )+alsa2y
        bbm(ny-1  )=bbm(ny-1  )+alsamy
        bbm(3     )=bbm(3     )+alsa3y
        bbm(ny-2  )=bbm(ny-2  )+alsaty
        bbm(4     )=bbm(4     )+alsa4y
        bbm(ny-3  )=bbm(ny-3  )+alsatty
        bbm(5:ny-4)=bbm(5:ny-4)+alsajy
     else
        bbm(2     )=bbm(2     )+alsa2y/pp2y(3)
        bbm(ny-1  )=bbm(ny-1  )+alsamy/pp2y(ny)
        bbm(3     )=bbm(3     )+alsa3y/pp2y(4)
        bbm(ny-2  )=bbm(ny-2  )+alsaty/pp2y(ny-1)
        bbm(4     )=bbm(4     )+alsa4y/pp2y(5)
        bbm(ny-3  )=bbm(ny-3  )+alsatty/pp2y(ny-2)
        bbm(5:ny-4)=bbm(5:ny-4)+alsajy/pp2y(6:ny-3)
     endif
     !BC on bbm
     bbm(1 )=zero
     bbm(ny)=zero
     !
     !DIAG SUP 2
     ccm(1     )=cs1y
     ccm(ny    )=csny
     ccm(2     )=0.!bs2y
     ccm(ny-1  )=0.!bsmy
     ccm(3     )=bs3y
     ccm(ny-2  )=bsty
     ccm(4     )=bs4y
     ccm(ny-3  )=bstty
     ccm(5:ny-4)=bsjy
     ccm = -xcst*ccm
     !BC on ccm
     ccm(1 )=zero
     ccm(ny)=zero
     !
     !DIAG SUP 3
     rrm(1     )=ds1y
     rrm(ny    )=dsny
     rrm(2     )=zero!cs2y
     rrm(ny-1  )=zero!csmy
     rrm(3     )=zero!cs3y
     rrm(ny-2  )=zero!csty
     rrm(4     )=cs4y
     rrm(ny-3  )=cstty
     rrm(5:ny-4)=csjy
     rrm = -xcst*rrm
     !BC on rrm
     rrm(1 )=zero
     rrm(ny)=zero
     !
     !DIAG SUP 4
     ttm(1     )=zero
     ttm(ny    )=zero
     ttm(2     )=zero
     ttm(ny-1  )=zero
     ttm(3     )=zero!ds3y
     ttm(ny-2  )=zero!dsty
     ttm(4     )=zero!ds4y
     ttm(ny-3  )=zero!dstty
     ttm(5:ny-4)=dsjy
     ttm = -xcst*ttm
     !BC on ttm
     ttm(1 )=zero
     ttm(ny)=zero
     !
     !DIAG INF 1
     if (istret==0) then
        ddm=bbm
     else
        ddm(1     )=bs1y
        ddm(ny    )=bsny
        ddm(2     )=as2y
        ddm(ny-1  )=asmy
        ddm(3     )=as3y
        ddm(ny-2  )=asty
        ddm(4     )=as4y
        ddm(ny-3  )=astty
        ddm(5:ny-4)=asjy
        ddm = -xcst*ddm
        ddm(2     )=ddm(2     )+alsa2y/pp2y(1)
        ddm(ny-1  )=ddm(ny-1  )+alsamy/pp2y(ny-2)
        ddm(3     )=ddm(3     )+alsa3y/pp2y(2)
        ddm(ny-2  )=ddm(ny-2  )+alsaty/pp2y(ny-3)
        ddm(4     )=ddm(4     )+alsa4y/pp2y(3)
        ddm(ny-3  )=ddm(ny-3  )+alsatty/pp2y(ny-4)
        ddm(5:ny-4)=ddm(5:ny-4)+alsajy/pp2y(4:ny-5)
        !BC on ddm
        ddm(1 )=zero
        ddm(ny)=zero
     endif
     !
     !DIAG INF 2
     eem=ccm
     !
     !DIAG INF 3
     qqm=rrm

     !DIAG INF 4
     uum=ttm

!!! NCL = 1, npaire=0, neumann BC, odd function
     !
     ! DIAG
     aam10(1     )=zero
     aam10(ny    )=zero
     aam10(2     )=-two*asjy-three*bsjy-two*csjy
     aam10(ny-1  )=aam10(2)
     aam10(3:ny-2)=-two*(asjy+bsjy+csjy)
     if (istret==0) then
        aam10 = one - xcst*aam10
     else
        aam10 = one/pp2y - xcst*aam10
     endif
     !
     !DIAG SUP 1
     bbm10(1     )=zero
     bbm10(ny    )=zero
     bbm10(2     )=asjy-csjy
     bbm10(ny-1  )=asjy
     bbm10(3     )=asjy
     bbm10(ny-2  )=asjy-csjy
     bbm10(4:ny-3)=asjy
     if (istret==0) then
        bbm10 = alsajy - xcst*bbm10
     else
        bbm10(2:ny-1) = alsajy/pp2y(3:ny) - xcst*bbm10(2:ny-1)
     endif
     !Bc on bbm10
     bbm10(1 )=zero
     bbm10(ny)=zero
     !
     !DIAG SUP 2
     ccm10(1     )=zero
     ccm10(ny    )=zero
     ccm10(2     )=bsjy
     ccm10(ny-1  )=zero
     ccm10(3     )=bsjy
     ccm10(ny-2  )=bsjy
     ccm10(4:ny-3)=bsjy
     ccm10 = -xcst*ccm10
     !
     !DIAG SUP 3
     rrm10(1     )=zero
     rrm10(ny    )=zero
     rrm10(2     )=csjy
     rrm10(ny-1  )=zero
     rrm10(3     )=csjy
     rrm10(ny-2  )=zero
     rrm10(4:ny-3)=csjy
     rrm10 = -xcst*rrm10
     !
     !DIAG INF 1
     ddm10(1     )=zero
     ddm10(ny    )=zero
     ddm10(2     )=asjy
     ddm10(ny-1  )=asjy-csjy
     ddm10(3     )=asjy-csjy
     ddm10(ny-2  )=asjy
     ddm10(4:ny-3)=asjy
     if (istret==0) then
        ddm10 = alsajy - xcst*ddm10
     else
        ddm10(2:ny-1) = alsajy/pp2y(1:ny-2) - xcst*ddm10(2:ny-1)
     endif
     !BC on ddm10
     ddm10(1 )=zero
     ddm10(ny)=zero
     !
     !DIAG INF 2
     eem10(1     )=zero
     eem10(ny    )=zero
     eem10(2     )=zero
     eem10(ny-1  )=bsjy
     eem10(3     )=bsjy
     eem10(ny-2  )=bsjy
     eem10(4:ny-3)=bsjy
     eem10 = -xcst*eem10
     !
     !DIAG INF 3
     qqm10(1     )=zero
     qqm10(ny    )=zero
     qqm10(2     )=zero
     qqm10(ny-1  )=csjy
     qqm10(3     )=zero
     qqm10(ny-2  )=csjy
     qqm10(4:ny-3)=csjy
     qqm10 = -xcst*qqm10

!!! NCL = 1, npaire=1, neumann BC, even function
     !
     ! DIAG
     aam11(1     )=-two*(asjy+bsjy+csjy)
     aam11(ny    )=aam11(1)
     aam11(2     )=-two*asjy-bsjy-two*csjy
     aam11(ny-1  )=aam11(2)
     aam11(3:ny-2)=-two*(asjy+bsjy+csjy)
     if (istret==0) then
        aam11 = one - xcst*aam11
     else
        aam11 = one/pp2y - xcst*aam11
     endif
     !
     !DIAG SUP 1
     bbm11(1     )=two*asjy
     bbm11(ny    )=zero
     bbm11(2     )=asjy+csjy
     bbm11(ny-1  )=asjy
     bbm11(3     )=asjy
     bbm11(ny-2  )=asjy+csjy
     bbm11(4:ny-3)=asjy
     if (istret==0) then
        bbm11 = alsajy - xcst*bbm11
     else
        bbm11(1:ny-1) = alsajy/pp2y(2:ny) - xcst*bbm11(1:ny-1)
     endif
     !BC on bbm11
     if (istret==0) then
        bbm11(1 )=bbm11(1)+alsajy
     else
        bbm11(1 )=bbm11(1)+alsajy/pp2y(2)
     endif
     bbm11(ny)=zero
     !
     !DIAG SUP 2
     ccm11(1     )=two*bsjy
     ccm11(ny    )=zero
     ccm11(2     )=bsjy
     ccm11(ny-1  )=zero
     ccm11(3     )=bsjy
     ccm11(ny-2  )=bsjy
     ccm11(4:ny-3)=bsjy
     ccm11 = -xcst*ccm11
     !
     !DIAG SUP 3
     rrm11(1     )=two*csjy
     rrm11(ny    )=zero
     rrm11(2     )=csjy
     rrm11(ny-1  )=zero
     rrm11(3     )=csjy
     rrm11(ny-2  )=zero
     rrm11(4:ny-3)=csjy
     rrm11 = -xcst*rrm11
     !
     !DIAG INF 1
     ddm11(1     )=zero
     ddm11(ny    )=two*asjy
     ddm11(2     )=asjy
     ddm11(ny-1  )=asjy+csjy
     ddm11(3     )=asjy+csjy
     ddm11(ny-2  )=asjy
     ddm11(4:ny-3)=asjy
     if (istret==0) then
        ddm11 = alsajy - xcst*ddm11
     else
        ddm11(2:ny) = alsajy/pp2y(1:ny-1) - xcst*ddm11(2:ny)
     endif
     !BC on ddm11
     ddm11(1 )=zero
     if (istret==0) then
        ddm11(ny)=ddm11(ny)+alsajy!a1
     else
        ddm11(ny)=ddm11(ny)+alsajy/pp2y(ny-1)!a1
     endif
     !
     !DIAG INF 2
     eem11(1     )=zero
     eem11(ny    )=two*bsjy
     eem11(2     )=zero
     eem11(ny-1  )=bsjy
     eem11(3     )=bsjy
     eem11(ny-2  )=bsjy
     eem11(4:ny-3)=bsjy
     eem11 = -xcst*eem11
     !
     !DIAG INF 3
     qqm11(1     )=zero
     qqm11(ny    )=two*csjy
     qqm11(2     )=zero
     qqm11(ny-1  )=csjy
     qqm11(3     )=zero
     qqm11(ny-2  )=csjy
     qqm11(4:ny-3)=csjy
     qqm11 = -xcst*qqm11

!!! NXL = 0
     !DIAG
     if (istret==0) then
        aam0 = one-xcst*(-two*(asjy+bsjy+csjy))
     else
        aam0 = one/pp2y-xcst*(-two*(asjy+bsjy+csjy))
     endif
     !
     !DIAG SUP 1
     if (istret==0) then
        bbm0 = alsajy-xcst*asjy
     else
        bbm0(1:ny-1) = alsajy/pp2y(2:ny) -xcst*asjy
        bbm0(ny) = alsajy/pp2y(1) -xcst*asjy
     endif
     !
     !DIAG SUP 2
     ccm0 = -xcst*bsjy
     !
     !DIAG SUP 3
     rrm0 = -xcst*csjy
     !
     !DIAG INF 1
     if (istret==0) then
        ddm0=bbm0
     else
        ddm0(1) = alsajy/pp2y(ny) -xcst*asjy
        ddm0(2:ny) = alsajy/pp2y(1:ny-1) -xcst*asjy
     endif
     !
     !DIAG INF 2
     eem0=ccm0
     !
     !DIAG INF 3
     qqm0=rrm0
  endif

  if (isecondder.ne.5) then
     call ludecomp7(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,&
          vvm,wwm,zzm,ny)
     call ludecomp7(aam10,bbm10,ccm10,ddm10,eem10,qqm10,ggm10,hhm10,ssm10,rrm10,&
          vvm10,wwm10,zzm10,ny)
     call ludecomp7(aam11,bbm11,ccm11,ddm11,eem11,qqm11,ggm11,hhm11,ssm11,rrm11,&
          vvm11,wwm11,zzm11,ny)
     call ludecomp7(aam0,bbm0,ccm0,ddm0,eem0,qqm0,ggm0,hhm0,ssm0,rrm0,&
          vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ny)
  else
     call ludecomp9(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ttm,uum,sssm,zzzm,ny)
     !NEED TO BE DONE: deal with other cases
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! END MATRIX M2 TIME IMPLICIT    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine implicit_schemes

!************************************************************
!
subroutine scalarimp(ux1,uy1,uz1,phi1,dphi1,is)
  !
  !************************************************************
  !IMPLICIT TIME INTEGRATION FOR D2/DY2
  !
  !************************************************************
  USE param
  USE variables
  USE var, ONLY: di1,tg1,th1,ti1,td1,ta2,tb2,tc2,td2,di2,di3,ta3,tb3
  USE decomp_2d
  USE derivY
  USE MPI
  use matinv

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphi1
  integer :: is

  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3
  real(mytype),dimension(ysize(1),ysize(3)) :: bc1
  
  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz,code

  !parametres CL
  real(mytype) :: x,y,z,r,lambda,phislbda,adiab,tjet,liss


  tg1=zero;th1=zero;ti1=zero;td1=zero
  ta2=zero;tb2=zero;tc2=zero;td2=zero
  ta3=zero;tb3=zero

  nvect1=xsize(1)*xsize(2)*xsize(3)
  nvect2=ysize(1)*ysize(2)*ysize(3)
  nvect3=zsize(1)*zsize(2)*zsize(3)
 
  !X PENCILS
  tg1=ux1*phi1
  call derx (th1,tg1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! npaire=0 pour phi
  call derxxS (tg1,phi1,di1,sx,sfxS,ssxS,swxS,xsize(1),xsize(2),xsize(3),0) ! npaire=0
  call transpose_x_to_y(phi1,phi2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(uz1,uz2)

  
  !Y PENCILS
  tc2=phi2*uy2
  call dery (tb2,tc2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=0 pour phi

  if (istret.ne.0) then

     call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0

     do k=1,ysize(3)
        do j=1,ysize(2)
           do i=1,ysize(1)
              ta2(i,j,k)=-pp4y(j)*tc2(i,j,k)
           enddo
        enddo
     enddo

  endif

  call transpose_y_to_z(phi2,phi3)
  call transpose_y_to_z(uz2,uz3)


  
  !Z PENCILS
  ta3=uz3*phi3
  call derz (tb3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! npaire=0 pour phi

  call derzzS (ta3,phi3,di3,sz,sfzS,sszS,swzS,zsize(1),zsize(2),zsize(3),0) ! npaire=0

  call transpose_z_to_y(ta3,tc2)
  call transpose_z_to_y(tb3,td2)

  !Y PENCILS ADD TERMS
  do k=1,ysize(3)
     do j=1,ysize(2)
        do i=1,ysize(1)
           tc2(i,j,k)=tc2(i,j,k)+ta2(i,j,k) !SECOND DERIVATIVE
           !         td2(i,j,k)=uz2(i,j,k)*td2(i,j,k)+uy2(i,j,k)*tb2(i,j,k) !FIRST DERIVATIVE
           td2(i,j,k)=td2(i,j,k)+tb2(i,j,k) !FIRST DERIVATIVE
        enddo
     enddo
  enddo

  call transpose_y_to_x(tc2,ti1)
  call transpose_y_to_x(td2,td1)


  !X PENCILS ADD TERMS
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           tg1(i,j,k)= tg1(i,j,k)+ti1(i,j,k) !SECOND DERIVATIVE
           !         th1(i,j,k)= ux1(i,j,k)*th1(i,j,k)+td1(i,j,k) !FIRST DERIVATIVE
           th1(i,j,k)= th1(i,j,k)+td1(i,j,k) !FIRST DERIVATIVE
        enddo
     enddo
  enddo

  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
           tg1(i,j,k)= (xnu/sc(is))*tg1(i,j,k)-th1(i,j,k)
        enddo
     enddo
  enddo

  ! TIME ADVANCEMENT EXPLICIT AB + IMPLICIT CN2 (d2/dy2)
  nxyz=xsize(1)*xsize(2)*xsize(3)



!!!!!!!!!!!!!!!!!!!!
  !CN2+AB3         !!!
!!!!!!!!!!!!!!!!!!!!
  if (itimescheme==7) then
     if ((itime.eq.1).and.(irestart.eq.0)) then
        if (nrank==0) print *,'Scalar start with Euler',itime
        !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES

        do k=1,xsize(3)
           do j=1,xsize(2)
              do i=1,xsize(1)
                 td1(i,j,k)= dt*tg1(i,j,k)
                 dphi1(i,j,k,2)= tg1(i,j,k)
              enddo
           enddo
        enddo

     else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES
        if  ((itime.eq.2).and.(irestart.eq.0)) then
           if (nrank==0) print *,'then with AB2',itime


           do k=1,xsize(3)
              do j=1,xsize(2)
                 do i=1,xsize(1)
                    td1(i,j,k)= onepfive*dt*tg1(i,j,k)-half*dt*dphi1(i,j,k,2)
                    dphi1(i,j,k,3)=dphi1(i,j,k,2)
                    dphi1(i,j,k,2)= tg1(i,j,k)
                 enddo
              enddo
           enddo

        else !FINALLY EXPLICIT AB3 + IMPLICIT CN SCHEMES

           do k=1,xsize(3)
              do j=1,xsize(2)
                 do i=1,xsize(1)
                    td1(i,j,k)= adt(itr)*tg1(i,j,k)+bdt(itr)*dphi1(i,j,k,2)&
                         +cdt(itr)*dphi1(i,j,k,3)
                    dphi1(i,j,k,3)=dphi1(i,j,k,2)
                    dphi1(i,j,k,2)= tg1(i,j,k)
                 enddo
              enddo
           enddo

        endif
     endif
  endif

  !Y-PENCIL FOR MATRIX INVERSION
  call transpose_x_to_y(phi1,phi2)
  call transpose_x_to_y(td1,ta2)


  !BC FOR THE SCALAR
  if (itype.eq.itype_tbl) then
     bc1(:,:)=phi2(:,ny-1,:)
     !in order to mimick a Neumann BC at the top of the domain for the TBL
  endif
  
  !ta2: A.T_hat
  !td2:(A+xcstB).Tn
  if (isecondder.ne.5) then
     call multmatrix7(td2,ta2,phi2,1)
  else
     call multmatrix9(td2,ta2,phi2)
  endif
  !right hand side
  ta2(:,:,:) = ta2(:,:,:) + td2(:,:,:)
 

  if (nclyS1.eq.2) then
     !BC at the bottom for the SCALAR
     ta2(:,1       ,:)=one
!     ta2(:,1,:)=zero
  endif

  if (nclySn.eq.2) then
     !BC at the top for the SCALAR
!     ta2(:,ysize(2),:)=one
     ta2(:,ysize(2),:)=bc1(:,:)
  endif

  !Inversion linear system Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
  !Inversion of the linear system Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
  !if secondder=5, we need nona inversion
  !if isecondder is not 5, we need septa inversion

  !NEED TO BE DONE: CASE 1-2 and 2-1!!!!
  if (isecondder.ne.5) then
     phi2=0.;
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        call septinv(phi2,ta2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        call septinv(phi2,ta2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call septinv(phi2,ta2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
     endif
  endif
    if (isecondder.eq.5) then
       phi2=0.;
       !TO BE DONE: Different types of BC
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        !NOT READY
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        !NOT READY
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call nonainv(phi2,ta2,ggm,hhm,ssm,sssm,ttm,zzzm,zzm,wwm,vvm,ysize(1),ysize(2),ysize(3))
     endif
  endif

  
  call transpose_y_to_x(phi2,phi1)

end subroutine scalarimp



subroutine init_implicit ()

  USE decomp_2d
  USE param
  USE variables
  implicit none

  allocate(aam(ny),bbm(ny),ccm(ny),ddm(ny),eem(ny),ggm(ny),hhm(ny),wwm(ny),zzm(ny))
  allocate(rrm(ny),qqm(ny),vvm(ny),ssm(ny))
  allocate(sssm(ny),zzzm(ny),ttm(ny),uum(ny))
  allocate(aam10(ny),bbm10(ny),ccm10(ny),ddm10(ny),eem10(ny),ggm10(ny),hhm10(ny),wwm10(ny),zzm10(ny))
  allocate(rrm10(ny),qqm10(ny),vvm10(ny),ssm10(ny))
  allocate(aam11(ny),bbm11(ny),ccm11(ny),ddm11(ny),eem11(ny),ggm11(ny),hhm11(ny),wwm11(ny),zzm11(ny))
  allocate(rrm11(ny),qqm11(ny),vvm11(ny),ssm11(ny))
  allocate(aam0(ny),bbm0(ny),ccm0(ny),ddm0(ny),eem0(ny),ggm0(ny),hhm0(ny),wwm0(ny),zzm0(ny))
  allocate(rrm0(ny),qqm0(ny),vvm0(ny),ssm0(ny),l1m(ny),l2m(ny),l3m(ny),u1m(ny),u2m(ny),u3m(ny))
  
end subroutine init_implicit
