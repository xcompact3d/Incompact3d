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
  ! Interface for septadiag / nonadiag LU decomposition used in implicit mode
  ! septadiag is for when isecondder is not equal 5 ('classic' second order schemes)
  ! nonadiag is for when isecondder is equal to 5 ('diffusive' second order schemes)
  ! ludecompX_0 is called when ncly=0 (when the matrix is cyclic)
  ! ludecompX_12 is called when ncly=1 or 2
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
    
    integer :: i,j,k
    integer, intent(in) :: ny
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm

    ! a=diag, b=diag+1, c=diag+2, r=diag+3
    !         d=diag-1, e=diag-2, q=diag-3
    ggm=zero;hhm=zero;ssm=zero ! U part (diag, diag+1, diag+2, diag+3=rrm)
    vvm=zero;wwm=zero;zzm=zero ! L part (diag=1, diag-3, diag-2, diag-1)

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

    integer :: i,j,k
    integer, intent(in) :: ny
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm
    real(mytype),dimension(ny), intent(out) :: l1m,l2m,l3m
    real(mytype),dimension(ny), intent(out) :: u1m,u2m,u3m

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

    integer :: i,j,k
    integer, intent(in) :: ny
    real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm,ttm,uum
    real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm,sssm
    real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm,zzzm

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

    integer :: i,j,k
    integer, intent(in) :: ny
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
  ! Xinv_0 is called when ncly=0 (when the matrix is cyclic)
  ! Xinv_12 is called when ncly=1 or 2
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

    integer :: i,j,k
    integer, intent(in) :: nx,ny,nz
    real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
    real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
    real(mytype),dimension(ny),intent(in)    :: ggm,hhm,ssm,rrm
    real(mytype),dimension(ny),intent(in)    :: vvm,wwm,zzm
    
    do k=1,nz
       !going down
       do i=1,nx
          xsol(i,1,k)=bbb(i,1,k)
       enddo
       do i=1,nx
          xsol(i,2,k)=bbb(i,2,k)-zzm(2)*xsol(i,1,k)
       enddo
       do i=1,nx
          xsol(i,3,k)=bbb(i,3,k)-zzm(3)*xsol(i,2,k)-wwm(3)*xsol(i,1,k)
       enddo

       do j=4,ny
          do i=1,nx
             xsol(i,j,k)=bbb(i,j,k)-vvm(j)*xsol(i,j-3,k)-wwm(j)*xsol(i,j-2,k) &
                  -zzm(j)*xsol(i,j-1,k)
          enddo
       enddo

       !going up
       do i=nx,1,-1
          xsol(i,ny,k)=xsol(i,ny,k)/ggm(ny)
       enddo
       do i=nx,1,-1
          xsol(i,ny-1,k)=(xsol(i,ny-1,k)-hhm(ny-1)*xsol(i,ny,k))/ggm(ny-1)
       enddo
       do i=nx,1,-1
          xsol(i,ny-2,k)=(xsol(i,ny-2,k)-hhm(ny-2)*xsol(i,ny-1,k)- &
               ssm(ny-2)*xsol(i,ny,k))/ggm(ny-2)
       enddo

       do j=ny-3,1,-1
          do i=nx,1,-1
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

    integer :: i,j,k,kk
    integer, intent(in) :: nx,ny,nz
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

    integer :: i,j,k
    integer, intent(in) :: nx,ny,nz
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

    integer :: i,j,k,kk
    integer, intent(in) :: nx,ny,nz
    integer :: code,ierror
    real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
    real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
    real(mytype),dimension(ny), intent(in)    :: ggm,hhm,ssm,rrm
    real(mytype),dimension(ny), intent(in)    :: vvm,wwm,zzm
    real(mytype),dimension(ny), intent(in) :: l1m,l2m,l3m
    real(mytype),dimension(ny), intent(in) :: u1m,u2m,u3m

    write(*,*) 'NOT READY YET! SIMULATION IS STOPPED!'
    call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop

  end subroutine nonainv_0
  
end module matinv

module ydiff_implicit

  implicit none

  private
  public :: inttimp, init_implicit, implicit_schemes

  contains
!
! Time integration, (semi)implicit Y diffusion
!    var1, input : variable at time n
!          output : variable at time n+1
!    dvar1 : r.h.s. of transport equation
!    npaire : odd / even variable, important when ncly*=1
!    isc : 0 for momentum, id of the scalar otherwise
!    forcing1 : r.h.s. term not present in dvar1 (pressure gradient)
!
subroutine  inttimp (var1,dvar1,npaire,isc,forcing1)

  USE MPI
  USE param
  USE variables
  USE var, ONLY: ta1, ta2, tb2, tc2, td2
  USE decomp_2d
  use derivY
  use matinv

  implicit none

  integer :: i,j,k,code,ierror

  !! IN
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)), optional, intent(in) :: forcing1
  integer, intent(in) :: npaire, isc

  !! IN/OUT
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

  !! LOCAL
  real(mytype),dimension(ysize(1),ysize(3)) :: bctop, bcbot

  if (itimescheme.eq.1) then
     !>>> Explicit Euler
     ta1(:,:,:) = gdt(1) * dvar1(:,:,:,1)

  elseif (itimescheme.eq.2) then
     !>>> AB2
     if ((itime.eq.1).and.(irestart.eq.0)) then
        !>>> Start with explicit Euler
        ta1(:,:,:) = dt*dvar1(:,:,:,1)

     else
        !>>> Then AB2
        ta1(:,:,:) = adt(1)*dvar1(:,:,:,1) + bdt(1)*dvar1(:,:,:,2)

     endif

     !save (N+Lxz)(n-1)
     dvar1(:,:,:,2) = dvar1(:,:,:,1)


  else if (itimescheme.eq.3) then
     !>>> AB3
     if ((itime.eq.1).and.(irestart.eq.0)) then
        !>>> Start with explicit Euler
        ta1(:,:,:) = dt*dvar1(:,:,:,1)

     else if ((itime.eq.2).and.(irestart.eq.0)) then
        !>>> Then AB2
        ta1(:,:,:) = onepfive*dt*dvar1(:,:,:,1) - half*dt*dvar1(:,:,:,2)

     else
        !>>> Then AB3
        ta1(:,:,:) = adt(1)*dvar1(:,:,:,1) + bdt(1)*dvar1(:,:,:,2) &
                   + cdt(1)*dvar1(:,:,:,3)

     endif

     !save (N+Lxz)(n-1)
     dvar1(:,:,:,3)=dvar1(:,:,:,2)
     dvar1(:,:,:,2)=dvar1(:,:,:,1)

  elseif (itimescheme.eq.4) then
     !>>> AB4
     if (nrank.eq.0) then
        print *, "AB4 not implemented!"
        STOP
     endif

  else
     !>>> We should not be here
     if (nrank.eq.0) then
        print *, "Unrecognised implicit itimescheme: ", itimescheme
     endif
     call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop

  endif

  if (present(forcing1)) then
     if ( (irestart.eq.1).or.(itime.gt.1) ) then
        ta1(:,:,:) = ta1(:,:,:) - forcing1(:,:,:)
     endif
  endif

  !Y-PENCIL FOR MATRIX INVERSION
  call transpose_x_to_y(var1,tb2)

  call transpose_x_to_y(ta1,ta2)

  !
  ! Prepare boundary conditions
  ! Used only when ncly*=2
  !    velocity
  !       isc = 0
  !       Dirichlet BC
  !    scalars
  !       1 <= isc <= numscalar
  !       Flexible BC, alpha_sc T +- beta_sc dT/dy = g_sc
  !
  ! Specific cases first
  ! This is the location for exotic / nonhomogeneous boundary conditions
  !
  if (itype.eq.itype_tbl .and. isc.eq.0) then
     bcbot(:,:) = zero
     bctop(:,:) = tb2(:,ny-1,:)
     !in order to mimick a Neumann BC at the top of the domain for the TBL
  !
  ! Generic homogeneous cases after
  !
  else if (isc.ne.0) then
     bcbot(:,:) = g_sc(isc, 1)
     bctop(:,:) = g_sc(isc, 2)
  else
     bcbot(:,:) = zero
     bctop(:,:) = zero
  endif
 
  !ta2: A.uhat
  !td2:(A+xcstB).un
  !if isecondder=5, we need nona inversion
  !id isecondder is not 5, we need septa inversion

  if (isecondder.ne.5) then
     if (isc.eq.0) then
        call multmatrix7(td2,ta2,tb2,npaire,ncly1,nclyn,xcst)
     else
        call multmatrix7(td2,ta2,tb2,npaire,nclyS1,nclySn,xcst_sc(isc))
     endif
  else if (isecondder.eq.5) then
     !TO BE DONE: Different types of BC
     if ((ncly1.eq.0).and.(nclyn.eq.0)) then
        !NOT READY
     elseif ((ncly1.eq.1).and.(nclyn.eq.1)) then
        !NOT READY
     elseif ((ncly1.eq.2).and.(nclyn.eq.2)) then
        call multmatrix9(td2,ta2,tb2,npaire)
     endif
  endif

  !Full second member BBB=A uhat+(A+xcst.B)u^n
  ta2(:,:,:)=td2(:,:,:)+ta2(:,:,:)

  !
  ! Apply boundary conditions
  !
  if ((isc.eq.0.and.ncly1.eq.2).or.(isc.gt.0.and.nclyS1.eq.2)) then
     ta2(:,1,:) = bcbot(:,:)
  endif
  if ((isc.eq.0.and.nclyn.eq.2).or.(isc.gt.0.and.nclySn.eq.2)) then
     ta2(:,ny,:) = bctop(:,:)
  endif
 
  !Inversion of the linear system Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
  !if secondder=5, we need nona inversion
  !if isecondder is not 5, we need septa inversion

  if (isecondder.ne.5) then
     if ((isc.eq.0).and.(ncly1.eq.0).and.(nclyn.eq.0)) then
        gg=>ggm0; hh=>hhm0; ss=>ssm0; rr=>rrm0; vv=>vvm0; ww=>wwm0; zz=>zzm0
        lo1=>l1m; lo2=>l2m; lo3=>l3m; up1=>u1m; up2=>u2m; up3=>u3m
     elseif ((isc.gt.0).and.(nclyS1.eq.0).and.(nclySn.eq.0)) then
        gg=>ggm0t(:,isc); hh=>hhm0t(:,isc); ss=>ssm0t(:,isc); rr=>rrm0t(:,isc); vv=>vvm0t(:,isc); ww=>wwm0t(:,isc); zz=>zzm0t(:,isc)
        lo1=>l1mt(:,isc); lo2=>l2mt(:,isc); lo3=>l3mt(:,isc); up1=>u1mt(:,isc); up2=>u2mt(:,isc); up3=>u3mt(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.1).and.(nclyn.eq.1).and.(npaire.eq.0)) then
        gg=>ggm10; hh=>hhm10; ss=>ssm10; rr=>rrm10; vv=>vvm10; ww=>wwm10; zz=>zzm10
     elseif ((isc.gt.0).and.(nclyS1.eq.1).and.(nclySn.eq.1).and.(npaire.eq.0)) then
        gg=>ggm10t(:,isc); hh=>hhm10t(:,isc); ss=>ssm10t(:,isc); rr=>rrm10t(:,isc); vv=>vvm10t(:,isc); ww=>wwm10t(:,isc); zz=>zzm10t(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.1).and.(nclyn.eq.1).and.(npaire.eq.1)) then
        gg=>ggm11; hh=>hhm11; ss=>ssm11; rr=>rrm11; vv=>vvm11; ww=>wwm11; zz=>zzm11
     elseif ((isc.gt.0).and.(nclyS1.eq.1).and.(nclySn.eq.1).and.(npaire.eq.1)) then
        gg=>ggm11t(:,isc); hh=>hhm11t(:,isc); ss=>ssm11t(:,isc); rr=>rrm11t(:,isc); vv=>vvm11t(:,isc); ww=>wwm11t(:,isc); zz=>zzm11t(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.2).and.(nclyn.eq.2)) then
        gg=>ggm; hh=>hhm; ss=>ssm; rr=>rrm; vv=>vvm; ww=>wwm; zz=>zzm
     elseif ((isc.gt.0).and.(nclyS1.eq.2).and.(nclySn.eq.2)) then
        gg=>ggmt(:,isc); hh=>hhmt(:,isc); ss=>ssmt(:,isc); rr=>rrmt(:,isc); vv=>vvmt(:,isc); ww=>wwmt(:,isc); zz=>zzmt(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.1).and.(nclyn.eq.2).and.(npaire.eq.0)) then
        gg=>ggm120; hh=>hhm120; ss=>ssm120; rr=>rrm120; vv=>vvm120; ww=>wwm120; zz=>zzm120
     elseif ((isc.gt.0).and.(nclyS1.eq.1).and.(nclySn.eq.2).and.(npaire.eq.0)) then
        gg=>ggm120t(:,isc); hh=>hhm120t(:,isc); ss=>ssm120t(:,isc); rr=>rrm120t(:,isc); vv=>vvm120t(:,isc); ww=>wwm120t(:,isc); zz=>zzm120t(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.1).and.(nclyn.eq.2).and.(npaire.eq.1)) then
        gg=>ggm121; hh=>hhm121; ss=>ssm121; rr=>rrm121; vv=>vvm121; ww=>wwm121; zz=>zzm121
     elseif ((isc.gt.0).and.(nclyS1.eq.1).and.(nclySn.eq.2).and.(npaire.eq.1)) then
        gg=>ggm121t(:,isc); hh=>hhm121t(:,isc); ss=>ssm121t(:,isc); rr=>rrm121t(:,isc); vv=>vvm121t(:,isc); ww=>wwm121t(:,isc); zz=>zzm121t(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.2).and.(nclyn.eq.1).and.(npaire.eq.0)) then
        gg=>ggm210; hh=>hhm210; ss=>ssm210; rr=>rrm210; vv=>vvm210; ww=>wwm210; zz=>zzm210
     elseif ((isc.gt.0).and.(nclyS1.eq.2).and.(nclySn.eq.1).and.(npaire.eq.0)) then
        gg=>ggm210t(:,isc); hh=>hhm210t(:,isc); ss=>ssm210t(:,isc); rr=>rrm210t(:,isc); vv=>vvm210t(:,isc); ww=>wwm210t(:,isc); zz=>zzm210t(:,isc)
     elseif ((isc.eq.0).and.(ncly1.eq.2).and.(nclyn.eq.1).and.(npaire.eq.1)) then
        gg=>ggm211; hh=>hhm211; ss=>ssm211; rr=>rrm211; vv=>vvm211; ww=>wwm211; zz=>zzm211
     elseif ((isc.gt.0).and.(nclyS1.eq.2).and.(nclySn.eq.1).and.(npaire.eq.1)) then
        gg=>ggm211t(:,isc); hh=>hhm211t(:,isc); ss=>ssm211t(:,isc); rr=>rrm211t(:,isc); vv=>vvm211t(:,isc); ww=>wwm211t(:,isc); zz=>zzm211t(:,isc)
     else
        ! We should not be here
        if (nrank == 0) then
           write(*,*)  "Error for time-implicit Y diffusion."
           if (isc == 0) write(*,*)  "   Wrong combination for ncly1, nclyn and npaire", ncly1, nclyn, npaire
           if (isc /= 0) write(*,*)  "   Wrong combination for nclyS1, nclySn and npaire", nclyS1, nclySn, npaire
        endif
        call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
     endif
     tb2=0.;
     if ((isc.eq.0.and.ncly1.eq.0.and.nclyn.eq.0).or.(isc.gt.0.and.nclyS1.eq.0.and.nclySn.eq.0)) then
        call septinv(tb2,ta2,gg,hh,ss,rr,vv,ww,zz,lo1,lo2,lo3,up1,up2,up3,ysize(1),ysize(2),ysize(3))
        nullify(lo1,lo2,lo3,up1,up2,up3)
     else
        call septinv(tb2,ta2,gg,hh,ss,rr,vv,ww,zz,ysize(1),ysize(2),ysize(3))
     endif
     nullify(gg,hh,ss,rr,vv,ww,zz)
  else if (isecondder.eq.5) then
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
  
  call transpose_y_to_x(tb2,var1)

  if (present(forcing1)) then
     if ( (irestart.eq.1).or.(itime.gt.1) ) then
        var1(:,:,:)=var1(:,:,:)+forcing1(:,:,:)
     endif
  endif


  return
end subroutine inttimp

!********************************************************************
!PREMULT FOR INTTIMP WITH SEPTADIAG
subroutine multmatrix7(td2,ta2,ux2,npaire,cly1,clyn,xcst)
! 
!********************************************************************
   USE param, only : iimplicit, istret, zero, ncly1, nclyn, two
   USE variables
   USE derivY
   USE decomp_2d
   USE ibm_param, only : ubcx,ubcy,ubcz
     
   implicit none
   
   integer, intent(in) :: npaire, cly1, clyn
   real(mytype), intent(in) :: xcst
   integer :: i,j,k,code
   real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
   real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
   real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2
   
   ! Compute A.ta2, store it in ta2
   if (istret.ne.0) then
      do j=1,ysize(2)
         ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
      enddo
   endif

   ! j=1,2,3
   if (cly1.eq.0) then
      td2(:,1,:) = alsajy*ta2(:,2,:) + ta2(:,1,:) + alsajy*ta2(:,ysize(2),:)
      do j=2,3
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
   else if (cly1.eq.1) then
      if (npaire.eq.0) then
         td2(:,1,:) = ta2(:,1,:)
      else
         td2(:,1,:) = two*alsajy*ta2(:,2,:) + ta2(:,1,:)
      endif
      do j=2,3
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
   else
      td2(:,1,:) = zero
      td2(:,2,:) = alsa2y*ta2(:,1,:) + ta2(:,2,:) + alsa2y*ta2(:,3,:)
      td2(:,3,:) = alsa3y*ta2(:,2,:) + ta2(:,3,:) + alsa3y*ta2(:,4,:)
   endif

   do j=4,ysize(2)-3
      td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
   enddo

   ! j=ny-2,ny-1,ny
   if (clyn.eq.0) then
      do j=ysize(2)-2,ysize(2)-1
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      td2(:,ysize(2),:) = alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) + alsajy*ta2(:,1,:)
   elseif (clyn.eq.1) then
      do j=ysize(2)-2,ysize(2)-1
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      if (npaire.eq.0) then
         td2(:,ysize(2),:) = ta2(:,ysize(2),:)
      else
         td2(:,ysize(2),:) = two*alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:)
      endif
   else
      td2(:,ysize(2)-2,:) = alsaty*ta2(:,ysize(2)-3,:) + ta2(:,ysize(2)-2,:) + alsaty*ta2(:,ysize(2)-1,:)
      td2(:,ysize(2)-1,:) = alsamy*ta2(:,ysize(2)-2,:) + ta2(:,ysize(2)-1,:) + alsamy*ta2(:,ysize(2),:)
      td2(:,ysize(2),:) = zero
   endif
   ta2 = td2

! Compute (A+nu*dt.B).un
! nu*dt*B.un first, needed for CN, not needed for backward Euler
   if (iimplicit.eq.1) then

      td2(:,:,:) = zero

   else if ((cly1.eq.1.or.clyn.eq.1) .and. npaire.eq.0) then

      ! Check if we are solving momentum or scalars
      if (cly1.eq.ncly1 .and. clyn.eq.nclyn) then
         call deryy(td2,ux2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0, ubcx)
      else
         call deryyS(td2,ux2,di2,sy,sfyS,ssyS,swyS,ysize(1),ysize(2),ysize(3),0, ubcx)
      endif
   else
      ! Check if we are solving momentum or scalars
      if (cly1.eq.ncly1 .and. clyn.eq.nclyn) then
         call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1, ubcx)
      else
         call deryyS(td2,ux2,di2,sy,sfypS,ssypS,swypS,ysize(1),ysize(2),ysize(3),1, ubcx)
      endif

   endif
   td2(:,:,:) = xcst * td2(:,:,:)

   if (istret.ne.0) then
      do j=1,ysize(2)
         ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
      enddo
   endif

   ! j=1,2,3
   if (cly1.eq.0) then
      td2(:,1,:) = alsajy*ux2(:,2,:) + ux2(:,1,:) + alsajy*ux2(:,ysize(2),:) &
                 + td2(:,1,:)
      do j=2,3
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + td2(:,j,:)
      enddo
   else if (cly1.eq.1) then
      if (npaire.eq.0) then
         td2(:,1,:) = ux2(:,1,:) + td2(:,1,:)
      else
         td2(:,1,:) = two*alsajy*ux2(:,2,:) + ux2(:,1,:) &
                    + td2(:,1,:)
      endif
      do j=2,3
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + td2(:,j,:)
      enddo
   else
      td2(:,1,:) = zero
      td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
                 + td2(:,2,:)
      td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
                 + td2(:,3,:)
   endif

   do j=4,ysize(2)-3
      td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                 + td2(:,j,:)
   enddo

   ! j=ny-2,ny-1,ny
   if (clyn.eq.0) then
      do j=ysize(2)-2,ysize(2)-1
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + td2(:,j,:)
      enddo
      td2(:,ysize(2),:) = alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajy*ux2(:,1,:) &
                        + td2(:,ysize(2),:)
   elseif (clyn.eq.1) then
      do j=ysize(2)-2,ysize(2)-1
         td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                    + td2(:,j,:)
      enddo
      if (npaire.eq.0) then
         td2(:,ysize(2),:) = ux2(:,ysize(2),:) + td2(:,ysize(2),:)
      else
         td2(:,ysize(2),:) = two*alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) &
                           + td2(:,ysize(2),:)
      endif
   else
      td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
                          + td2(:,ysize(2)-2,:)
      td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
                          + td2(:,ysize(2)-1,:)
      td2(:,ysize(2),:) = zero
   endif

end subroutine multmatrix7

!********************************************************************

subroutine multmatrix9(td2,ta2,ux2,npaire)
  !
  !********************************************************************
  USE param
  USE variables
  USE derivY
  USE decomp_2d
  USE ibm_param, only : ubcx,ubcy,ubcz
  
  implicit none

  integer, intent(in) :: npaire
  integer :: i,j,k,code
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

     td2(:,1,:) = zero
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

  if (iimplicit.eq.1) then

     td2(:,:,:) = zero

  elseif ((ncly1.eq.1.or.nclyn.eq.1) .and. npaire.eq.0) then

     call deryy(td2,ux2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0, ubcx)
  else
     call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1, ubcx)         

  endif
  td2(:,:,:) = xcst * td2(:,:,:)

  if (istret.ne.0) then
     do j=1,ysize(2)
        ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
     enddo
  endif

  if (ncly1.eq.0) then

     td2(:,1,:) = alsajy*ux2(:,2,:) + ux2(:,1,:) + alsajy*ux2(:,ysize(2),:) &
                + td2(:,1,:)
     do j=2,ysize(2)-1
        td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                   + td2(:,j,:)
     enddo
     td2(:,ysize(2),:) = alsajy*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajy*ux2(:,1,:) &
                       + td2(:,ysize(2),:)

  elseif (ncly1.eq.1) then

  elseif (ncly1.eq.2) then

     td2(:,1,:) = zero
     td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
                + td2(:,2,:)
     td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
                + td2(:,3,:)
     td2(:,4,:) = alsa4y*ux2(:,3,:) + ux2(:,4,:) + alsa4y*ux2(:,5,:) &
                + td2(:,4,:)
     do j=5,ysize(2)-4
        td2(:,j,:) = alsajy*ux2(:,j-1,:) + ux2(:,j,:) + alsajy*ux2(:,j+1,:) &
                   + td2(:,j,:)
     enddo
     td2(:,ysize(2)-3,:) = alsatty*ux2(:,ysize(2)-4,:) + ux2(:,ysize(2)-3,:) + alsatty*ux2(:,ysize(2)-2,:) &
                         + td2(:,ysize(2)-3,:)
     td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
                         + td2(:,ysize(2)-2,:)
     td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
                         + td2(:,ysize(2)-1,:)
     td2(:,ysize(2),:) = zero

  endif

end subroutine multmatrix9

!
! Compute 1D arrays containing LU decomposition
!
subroutine implicit_schemes()

  USE param
  USE derivY
  USE variables
  USE var
  USE param
  use ludecomp

  implicit none

  integer  :: i,j,k,is

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MATRIX M2 TIME IMPLICIT !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!FOR VELOCITY  !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

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
     call init_implicit_coef(aam, aamt)
     if (istret==0) then
        aam = one-xcst*aam
        do is = 1, numscalar
           aamt(:,is) = one-xcst_sc(is)*aamt(:,is)
        enddo
     else
        aam = one/pp2y -xcst*aam
        do is = 1, numscalar
           aamt(:,is) = one/pp2y -xcst_sc(is)*aamt(:,is)
        enddo
     endif
     ! BC on aam, dirichlet
     aam(1 )=one
     aam(ny)=one
     ! BC on aamt
     if (istret.eq.0) then
        do is = 1, numscalar
           aamt(1 ,is) = alpha_sc(is,1) + beta_sc(is,1)*(11.d0/6.d0/dy)
           aamt(ny,is) = alpha_sc(is,2) + beta_sc(is,2)*(11.d0/6.d0/dy)
        enddo
     else
        do is = 1, numscalar
           aamt(1 ,is) = alpha_sc(is,1) + beta_sc(is,1)*ppy(1 )*(11.d0/6.d0/dy)
           aamt(ny,is) = alpha_sc(is,2) + beta_sc(is,2)*ppy(ny)*(11.d0/6.d0/dy)
        enddo
     endif
     !DIAG SUP 1
     bbm(1     )=bs1y
     bbm(ny    )=bsny
     bbm(2     )=as2y
     bbm(ny-1  )=asmy
     bbm(3     )=as3y
     bbm(ny-2  )=asty
     bbm(4:ny-3)=asjy
     call init_implicit_coef(bbm, bbmt)
     bbm = -xcst*bbm
     do is = 1, numscalar
        bbmt(:,is) = -xcst_sc(is)*bbmt(:,is)
     enddo
     if (istret==0) then
        bbm(2     )=bbm(2     )+alsa2y
        bbm(ny-1  )=bbm(ny-1  )+alsamy
        bbm(3     )=bbm(3     )+alsa3y
        bbm(ny-2  )=bbm(ny-2  )+alsaty
        bbm(4:ny-3)=bbm(4:ny-3)+alsajy
        do is = 1, numscalar
           bbmt(2     ,is)=bbmt(2     ,is)+alsa2y
           bbmt(ny-1  ,is)=bbmt(ny-1  ,is)+alsamy
           bbmt(3     ,is)=bbmt(3     ,is)+alsa3y
           bbmt(ny-2  ,is)=bbmt(ny-2  ,is)+alsaty
           bbmt(4:ny-3,is)=bbmt(4:ny-3,is)+alsajy
        enddo
     else
        bbm(2     )=bbm(2     )+alsa2y/pp2y(3)
        bbm(ny-1  )=bbm(ny-1  )+alsamy/pp2y(ny)
        bbm(3     )=bbm(3     )+alsa3y/pp2y(4)
        bbm(ny-2  )=bbm(ny-2  )+alsaty/pp2y(ny-1)
        bbm(4:ny-3)=bbm(4:ny-3)+alsajy/pp2y(5:ny-2)
        do is = 1, numscalar
           bbmt(2     ,is)=bbmt(2     ,is)+alsa2y/pp2y(3)
           bbmt(ny-1  ,is)=bbmt(ny-1  ,is)+alsamy/pp2y(ny)
           bbmt(3     ,is)=bbmt(3     ,is)+alsa3y/pp2y(4)
           bbmt(ny-2  ,is)=bbmt(ny-2  ,is)+alsaty/pp2y(ny-1)
           bbmt(4:ny-3,is)=bbmt(4:ny-3,is)+alsajy/pp2y(5:ny-2)
        enddo
     endif
     ! BC on bbm, dirichlet
     bbm(1 )=zero
     bbm(ny)=zero
     ! BC on bbmt
     if (istret.eq.0) then
        do is = 1, numscalar
           bbmt(1 ,is) = beta_sc(is,1)*(-18.d0/6.d0/dy)
           bbmt(ny,is) = zero
        enddo
     else
        do is = 1, numscalar
           bbmt(1 ,is) = beta_sc(is,1)*ppy(1)*(-18.d0/6.d0/dy)
           bbmt(ny,is) = zero
        enddo
     endif
     !DIAG SUP 2
     ccm(1     )=cs1y
     ccm(ny    )=csny
     ccm(2     )=zero
     ccm(ny-1  )=zero
     ccm(3     )=bs3y
     ccm(ny-2  )=bsty
     ccm(4:ny-3)=bsjy
     call init_implicit_coef(ccm, ccmt)
     ccm = -xcst*ccm
     do is = 1, numscalar
        ccmt(:,is) = -xcst_sc(is)*ccmt(:,is)
     enddo
     ! BC on ccm
     ccm(1 )=zero
     ccm(ny)=zero
     ! BC on ccmt
     if (istret.eq.0) then
        do is = 1, numscalar
           ccmt(1 ,is) = beta_sc(is,1)*(nine/six/dy)
           ccmt(ny,is) = zero
        enddo
     else
        do is = 1, numscalar
           ccmt(1 ,is) = beta_sc(is,1)*ppy(1)*(nine/six/dy)
           ccmt(ny,is) = zero
        enddo
     endif
     !DIAG SUP 3
     rrm(1     )=ds1y
     rrm(ny    )=dsny
     rrm(2     )=zero
     rrm(ny-1  )=zero
     rrm(3     )=zero
     rrm(ny-2  )=zero
     rrm(4:ny-3)=csjy
     call init_implicit_coef(rrm, rrmt)
     rrm = -xcst*rrm
     do is = 1, numscalar
        rrmt(:,is) = -xcst_sc(is)*rrmt(:,is)
     enddo
     ! BC on rrm
     rrm(1 )=zero
     rrm(ny)=zero
     ! BC on rrmt
     if (istret.eq.0) then
        do is = 1, numscalar
           rrmt(1 ,is) = beta_sc(is,1)*(-two/six/dy)
           rrmt(ny,is) = zero
        enddo
     else
        do is = 1, numscalar
           rrmt(1 ,is) = beta_sc(is,1)*ppy(1)*(-two/six/dy)
           rrmt(ny,is) = zero
        enddo
     endif
     !DIAG INF 1
     if (istret==0) then
        ddm=bbm
        ddmt = bbmt
     else
        ddm(1     )=bs1y
        ddm(ny    )=bsny
        ddm(2     )=as2y
        ddm(ny-1  )=asmy
        ddm(3     )=as3y
        ddm(ny-2  )=asty
        ddm(4:ny-3)=asjy
        call init_implicit_coef(ddm, ddmt)
        ddm = -xcst*ddm
        do is = 1, numscalar
           ddmt(:,is) = -xcst_sc(is)*ddmt(:,is)
        enddo
        ddm(2     )=ddm(2     )+alsa2y/pp2y(1)
        ddm(ny-1  )=ddm(ny-1  )+alsamy/pp2y(ny-2)
        ddm(3     )=ddm(3     )+alsa3y/pp2y(2)
        ddm(ny-2  )=ddm(ny-2  )+alsaty/pp2y(ny-3)
        ddm(4:ny-3)=ddm(4:ny-3)+alsajy/pp2y(3:ny-4)
        do is = 1, numscalar
           ddmt(2     ,is)=ddmt(2     ,is)+alsa2y/pp2y(1)
           ddmt(ny-1  ,is)=ddmt(ny-1  ,is)+alsamy/pp2y(ny-2)
           ddmt(3     ,is)=ddmt(3     ,is)+alsa3y/pp2y(2)
           ddmt(ny-2  ,is)=ddmt(ny-2  ,is)+alsaty/pp2y(ny-3)
           ddmt(4:ny-3,is)=ddmt(4:ny-3,is)+alsajy/pp2y(3:ny-4)
        enddo
        ! BC on ddm
        ddm(1 )=zero
        ddm(ny)=zero
     endif
     ! BC on ddmt
     if (istret.eq.0) then
        do is = 1, numscalar
           ddmt(1 ,is) = zero
           ddmt(ny,is) = beta_sc(is,2)*(-eighteen/six/dy)
        enddo
     else
        do is = 1, numscalar
           ddmt(1 ,is) = zero
           ddmt(ny,is) = beta_sc(is,2)*ppy(ny)*(-eighteen/six/dy)
        enddo
     endif
     !
     !DIAG INF 2
     eem=ccm
     eemt = ccmt
     ! BC on eemt
     if (istret.eq.0) then
        do is = 1, numscalar
           eemt(1 ,is) = zero
           eemt(ny,is) = beta_sc(is,2)*(nine/six/dy)
        enddo
     else
        do is = 1, numscalar
           eemt(1 ,is) = zero
           eemt(ny,is) = beta_sc(is,2)*ppy(ny)*(nine/six/dy)
        enddo
     endif
     !
     !DIAG INF 3
     qqm=rrm
     qqmt = rrmt
     ! BC on rrmt
     if (istret.eq.0) then
        do is = 1, numscalar
           qqmt(1 ,is) = zero
           qqmt(ny,is) = beta_sc(is,2)*(-two/six/dy)
        enddo
     else
        do is = 1, numscalar
           qqmt(1 ,is) = zero
           qqmt(ny,is) = beta_sc(is,2)*ppy(ny)*(-two/six/dy)
        enddo
     endif

!!! NCL = 1, npaire=0, neumann BC, odd function
     !
     ! DIAG
     aam10(1     )=zero
     aam10(ny    )=zero
     aam10(2     )=-two*asjy-three*bsjy-two*csjy
     aam10(ny-1  )=aam10(2)
     aam10(3:ny-2)=-two*(asjy+bsjy+csjy)
     call init_implicit_coef(aam10, aam10t)
     if (istret==0) then
        aam10 = one - xcst*aam10
        do is = 1, numscalar
           aam10t(:,is) = one - xcst_sc(is)*aam10t(:,is)
        enddo
     else
        aam10 = one/pp2y - xcst*aam10
        do is = 1, numscalar
           aam10t(:,is) = one/pp2y - xcst_sc(is)*aam10t(:,is)
        enddo
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
     call init_implicit_coef(bbm10, bbm10t)
     if (istret==0) then
        bbm10 = alsajy - xcst*bbm10
        do is = 1, numscalar
           bbm10t(:,is) = alsajy - xcst_sc(is)*bbm10t(:,is)
        enddo
     else
        bbm10(2:ny-1) = alsajy/pp2y(3:ny) - xcst*bbm10(2:ny-1)
        do is = 1, numscalar
           bbm10t(2:ny-1,is) = alsajy/pp2y(3:ny) - xcst_sc(is)*bbm10t(2:ny-1,is)
        enddo
     endif
     !BC on bbm10
     bbm10(1 )=zero
     bbm10(ny)=zero
     ! BC on bbm10t
     bbm10t(1 ,:)=zero
     bbm10t(ny,:)=zero
     !
     !DIAG SUP 2
     ccm10(1     )=zero
     ccm10(ny    )=zero
     ccm10(2     )=bsjy
     ccm10(ny-1  )=zero
     ccm10(3     )=bsjy
     ccm10(ny-2  )=bsjy
     ccm10(4:ny-3)=bsjy
     call init_implicit_coef(ccm10, ccm10t)
     ccm10 = -xcst*ccm10
     do is = 1, numscalar
        ccm10t(:,is) = -xcst_sc(is)*ccm10t(:,is)
     enddo
     !
     !DIAG SUP 3
     rrm10(1     )=zero
     rrm10(ny    )=zero
     rrm10(2     )=csjy
     rrm10(ny-1  )=zero
     rrm10(3     )=csjy
     rrm10(ny-2  )=zero
     rrm10(4:ny-3)=csjy
     call init_implicit_coef(rrm10, rrm10t)
     rrm10 = -xcst*rrm10
     do is = 1, numscalar
        rrm10t(:,is) = -xcst_sc(is)*rrm10t(:,is)
     enddo
     !
     !DIAG INF 1
     ddm10(1     )=zero
     ddm10(ny    )=zero
     ddm10(2     )=asjy
     ddm10(ny-1  )=asjy-csjy
     ddm10(3     )=asjy-csjy
     ddm10(ny-2  )=asjy
     ddm10(4:ny-3)=asjy
     call init_implicit_coef(ddm10, ddm10t)
     if (istret==0) then
        ddm10 = alsajy - xcst*ddm10
        do is = 1, numscalar
           ddm10t(:,is) = alsajy - xcst_sc(is)*ddm10t(:,is)
        enddo
     else
        ddm10(2:ny-1) = alsajy/pp2y(1:ny-2) - xcst*ddm10(2:ny-1)
        do is = 1, numscalar
           ddm10t(2:ny-1,is) = alsajy/pp2y(1:ny-2) - xcst_sc(is)*ddm10t(2:ny-1,is)
        enddo
     endif
     !BC on ddm10
     ddm10(1 )=zero
     ddm10(ny)=zero
     ! BC on ddm10t
     ddm10t(1 ,:)=zero
     ddm10t(ny,:)=zero
     !
     !DIAG INF 2
     eem10(1     )=zero
     eem10(ny    )=zero
     eem10(2     )=zero
     eem10(ny-1  )=bsjy
     eem10(3     )=bsjy
     eem10(ny-2  )=bsjy
     eem10(4:ny-3)=bsjy
     call init_implicit_coef(eem10, eem10t)
     eem10 = -xcst*eem10
     do is = 1, numscalar
        eem10t(:,is) = -xcst_sc(is)*eem10t(:,is)
     enddo
     !
     !DIAG INF 3
     qqm10(1     )=zero
     qqm10(ny    )=zero
     qqm10(2     )=zero
     qqm10(ny-1  )=csjy
     qqm10(3     )=zero
     qqm10(ny-2  )=csjy
     qqm10(4:ny-3)=csjy
     call init_implicit_coef(qqm10, qqm10t)
     qqm10 = -xcst*qqm10
     do is = 1, numscalar
        qqm10t(:,is) = -xcst_sc(is)*qqm10t(:,is)
     enddo

!!! NCL = 1, npaire=1, neumann BC, even function
     !
     ! DIAG
     aam11(1     )=-two*(asjy+bsjy+csjy)
     aam11(ny    )=aam11(1)
     aam11(2     )=-two*asjy-bsjy-two*csjy
     aam11(ny-1  )=aam11(2)
     aam11(3:ny-2)=-two*(asjy+bsjy+csjy)
     call init_implicit_coef(aam11, aam11t)
     if (istret==0) then
        aam11 = one - xcst*aam11
        do is = 1, numscalar
           aam11t(:,is) = one - xcst_sc(is)*aam11t(:,is)
        enddo
     else
        aam11 = one/pp2y - xcst*aam11
        do is = 1, numscalar
           aam11t(:,is) = one/pp2y - xcst_sc(is)*aam11t(:,is)
        enddo
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
     call init_implicit_coef(bbm11, bbm11t)
     if (istret==0) then
        bbm11 = alsajy - xcst*bbm11
        do is = 1, numscalar
           bbm11t(:,is) = alsajy - xcst_sc(is)*bbm11t(:,is)
        enddo
     else
        bbm11(1:ny-1) = alsajy/pp2y(2:ny) - xcst*bbm11(1:ny-1)
        do is = 1, numscalar
           bbm11t(1:ny-1,is) = alsajy/pp2y(2:ny) - xcst_sc(is)*bbm11t(1:ny-1,is)
        enddo
     endif
     !BC on bbm11 and bbm11t
     if (istret==0) then
        bbm11(1 )=bbm11(1)+alsajy
        bbm11t(1 ,:)=bbm11t(1,:)+alsajy
     else
        bbm11(1 )=bbm11(1)+alsajy/pp2y(2)
        bbm11t(1 ,:)=bbm11t(1,:)+alsajy/pp2y(2)
     endif
     bbm11(ny)=zero
     bbm11t(ny,:)=zero
     !
     !DIAG SUP 2
     ccm11(1     )=two*bsjy
     ccm11(ny    )=zero
     ccm11(2     )=bsjy
     ccm11(ny-1  )=zero
     ccm11(3     )=bsjy
     ccm11(ny-2  )=bsjy
     ccm11(4:ny-3)=bsjy
     call init_implicit_coef(ccm11, ccm11t)
     ccm11 = -xcst*ccm11
     do is = 1, numscalar
        ccm11t(:,is) = -xcst_sc(is)*ccm11t(:,is)
     enddo
     !
     !DIAG SUP 3
     rrm11(1     )=two*csjy
     rrm11(ny    )=zero
     rrm11(2     )=csjy
     rrm11(ny-1  )=zero
     rrm11(3     )=csjy
     rrm11(ny-2  )=zero
     rrm11(4:ny-3)=csjy
     call init_implicit_coef(rrm11, rrm11t)
     rrm11 = -xcst*rrm11
     do is = 1, numscalar
        rrm11t(:,is) = -xcst_sc(is)*rrm11t(:,is)
     enddo
     !
     !DIAG INF 1
     ddm11(1     )=zero
     ddm11(ny    )=two*asjy
     ddm11(2     )=asjy
     ddm11(ny-1  )=asjy+csjy
     ddm11(3     )=asjy+csjy
     ddm11(ny-2  )=asjy
     ddm11(4:ny-3)=asjy
     call init_implicit_coef(ddm11, ddm11t)
     if (istret==0) then
        ddm11 = alsajy - xcst*ddm11
        do is = 1, numscalar
           ddm11t(:,is) = alsajy - xcst_sc(is)*ddm11t(:,is)
        enddo
     else
        ddm11(2:ny) = alsajy/pp2y(1:ny-1) - xcst*ddm11(2:ny)
        do is = 1, numscalar
           ddm11t(2:ny,is) = alsajy/pp2y(1:ny-1) - xcst_sc(is)*ddm11t(2:ny,is)
        enddo
     endif
     !BC on ddm11 and ddm11t
     ddm11(1 )=zero
     ddm11t(1 ,:)=zero
     if (istret==0) then
        ddm11(ny)=ddm11(ny)+alsajy!a1
        ddm11t(ny,:)=ddm11t(ny,:)+alsajy!a1
     else
        ddm11(ny)=ddm11(ny)+alsajy/pp2y(ny-1)!a1
        ddm11t(ny,:)=ddm11t(ny,:)+alsajy/pp2y(ny-1)!a1
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
     call init_implicit_coef(eem11, eem11t)
     eem11 = -xcst*eem11
     do is = 1, numscalar
        eem11t(:,is) = -xcst_sc(is)*eem11t(:,is)
     enddo
     !
     !DIAG INF 3
     qqm11(1     )=zero
     qqm11(ny    )=two*csjy
     qqm11(2     )=zero
     qqm11(ny-1  )=csjy
     qqm11(3     )=zero
     qqm11(ny-2  )=csjy
     qqm11(4:ny-3)=csjy
     call init_implicit_coef(qqm11, qqm11t)
     qqm11 = -xcst*qqm11
     do is = 1, numscalar
        qqm11t(:,is) = -xcst_sc(is)*qqm11t(:,is)
     enddo

!!! NXL = 0
     !DIAG
     if (istret==0) then
        aam0 = one-xcst*(-two*(asjy+bsjy+csjy))
        do is = 1, numscalar
           aam0t(:,is) = one-xcst_sc(is)*(-two*(asjy+bsjy+csjy))
        enddo
     else
        aam0 = one/pp2y-xcst*(-two*(asjy+bsjy+csjy))
        do is = 1, numscalar
           aam0t(:,is) = one/pp2y-xcst_sc(is)*(-two*(asjy+bsjy+csjy))
        enddo
     endif
     !
     !DIAG SUP 1
     if (istret==0) then
        bbm0 = alsajy-xcst*asjy
        do is = 1, numscalar
           bbm0t(:,is) = alsajy-xcst_sc(is)*asjy
        enddo
     else
        bbm0(1:ny-1) = alsajy/pp2y(2:ny) -xcst*asjy
        bbm0(ny) = alsajy/pp2y(1) -xcst*asjy
        do is = 1, numscalar
           bbm0t(1:ny-1,is) = alsajy/pp2y(2:ny) -xcst_sc(is)*asjy
           bbm0t(ny,is) = alsajy/pp2y(1) -xcst_sc(is)*asjy
        enddo
     endif
     !
     !DIAG SUP 2
     ccm0 = -xcst*bsjy
     do is = 1, numscalar
        ccm0t(:,is) = -xcst_sc(is)*bsjy
     enddo
     !
     !DIAG SUP 3
     rrm0 = -xcst*csjy
     do is = 1, numscalar
        rrm0t(:,is) = -xcst_sc(is)*csjy
     enddo
     !
     !DIAG INF 1
     if (istret==0) then
        ddm0=bbm0
        ddm0t = bbm0t
     else
        ddm0(1) = alsajy/pp2y(ny) -xcst*asjy
        ddm0(2:ny) = alsajy/pp2y(1:ny-1) -xcst*asjy
        do is = 1, numscalar
           ddm0t(1,is) = alsajy/pp2y(ny) -xcst_sc(is)*asjy
           ddm0t(2:ny,is) = alsajy/pp2y(1:ny-1) -xcst_sc(is)*asjy
        enddo
     endif
     !
     !DIAG INF 2
     eem0=ccm0
     eem0t=ccm0t
     !
     !DIAG INF 3
     qqm0=rrm0
     qqm0t=rrm0t

     ! velocity, ncly1 = 1, nclyn = 2, npaire = 0
     aam120=aam10; aam120(ny-2:ny)=aam(ny-2:ny)
     bbm120=bbm10; bbm120(ny-2:ny)=bbm(ny-2:ny)
     ccm120=ccm10; ccm120(ny-2:ny)=ccm(ny-2:ny)
     ddm120=ddm10; ddm120(ny-2:ny)=ddm(ny-2:ny)
     eem120=eem10; eem120(ny-2:ny)=eem(ny-2:ny)
     qqm120=qqm10; qqm120(ny-2:ny)=qqm(ny-2:ny)
     rrm120=rrm10; rrm120(ny-2:ny)=rrm(ny-2:ny)
     ! velocity, ncly1 = 1, nclyn = 2, npaire = 1
     aam121=aam11; aam121(ny-2:ny)=aam(ny-2:ny)
     bbm121=bbm11; bbm121(ny-2:ny)=bbm(ny-2:ny)
     ccm121=ccm11; ccm121(ny-2:ny)=ccm(ny-2:ny)
     ddm121=ddm11; ddm121(ny-2:ny)=ddm(ny-2:ny)
     eem121=eem11; eem121(ny-2:ny)=eem(ny-2:ny)
     qqm121=qqm11; qqm121(ny-2:ny)=qqm(ny-2:ny)
     rrm121=rrm11; rrm121(ny-2:ny)=rrm(ny-2:ny)
     ! velocity, ncly1 = 2, nclyn = 1, npaire = 0
     aam210=aam; aam210(ny-2:ny)=aam10(ny-2:ny)
     bbm210=bbm; bbm210(ny-2:ny)=bbm10(ny-2:ny)
     ccm210=ccm; ccm210(ny-2:ny)=vvm10(ny-2:ny)
     ddm210=ddm; ddm210(ny-2:ny)=ddm10(ny-2:ny)
     eem210=eem; eem210(ny-2:ny)=eem10(ny-2:ny)
     qqm210=qqm; qqm210(ny-2:ny)=qqm10(ny-2:ny)
     rrm210=rrm; rrm210(ny-2:ny)=rrm10(ny-2:ny)
     ! velocity, ncly1 = 2, nclyn = 1, npaire = 1
     aam211=aam; aam211(ny-2:ny)=aam11(ny-2:ny)
     bbm211=bbm; bbm211(ny-2:ny)=bbm11(ny-2:ny)
     ccm211=ccm; ccm211(ny-2:ny)=ccm11(ny-2:ny)
     ddm211=ddm; ddm211(ny-2:ny)=ddm11(ny-2:ny)
     eem211=eem; eem211(ny-2:ny)=eem11(ny-2:ny)
     qqm211=qqm; qqm211(ny-2:ny)=qqm11(ny-2:ny)
     rrm211=rrm; rrm211(ny-2:ny)=rrm11(ny-2:ny)

     ! scalars, ncly1 = 1, nclyn = 2, npaire = 0
     aam120t=aam10t; aam120t(ny-2:ny,:)=aamt(ny-2:ny,:)
     bbm120t=bbm10t; bbm120t(ny-2:ny,:)=bbmt(ny-2:ny,:)
     ccm120t=ccm10t; ccm120t(ny-2:ny,:)=ccmt(ny-2:ny,:)
     ddm120t=ddm10t; ddm120t(ny-2:ny,:)=ddmt(ny-2:ny,:)
     eem120t=eem10t; eem120t(ny-2:ny,:)=eemt(ny-2:ny,:)
     qqm120t=qqm10t; qqm120t(ny-2:ny,:)=qqmt(ny-2:ny,:)
     rrm120t=rrm10t; rrm120t(ny-2:ny,:)=rrmt(ny-2:ny,:)
     ! scalars, ncly1 = 1, nclyn = 2, npaire = 1
     aam121t=aam11t; aam121t(ny-2:ny,:)=aamt(ny-2:ny,:)
     bbm121t=bbm11t; bbm121t(ny-2:ny,:)=bbmt(ny-2:ny,:)
     ccm121t=ccm11t; ccm121t(ny-2:ny,:)=ccmt(ny-2:ny,:)
     ddm121t=ddm11t; ddm121t(ny-2:ny,:)=ddmt(ny-2:ny,:)
     eem121t=eem11t; eem121t(ny-2:ny,:)=eemt(ny-2:ny,:)
     qqm121t=qqm11t; qqm121t(ny-2:ny,:)=qqmt(ny-2:ny,:)
     rrm121t=rrm11t; rrm121t(ny-2:ny,:)=rrmt(ny-2:ny,:)
     ! scalars, ncly1 = 2, nclyn = 1, npaire = 0
     aam210t=aamt; aam210t(ny-2:ny,:)=aam10t(ny-2:ny,:)
     bbm210t=bbmt; bbm210t(ny-2:ny,:)=bbm10t(ny-2:ny,:)
     ccm210t=ccmt; ccm210t(ny-2:ny,:)=vvm10t(ny-2:ny,:)
     ddm210t=ddmt; ddm210t(ny-2:ny,:)=ddm10t(ny-2:ny,:)
     eem210t=eemt; eem210t(ny-2:ny,:)=eem10t(ny-2:ny,:)
     qqm210t=qqmt; qqm210t(ny-2:ny,:)=qqm10t(ny-2:ny,:)
     rrm210t=rrmt; rrm210t(ny-2:ny,:)=rrm10t(ny-2:ny,:)
     ! scalars, ncly1 = 2, nclyn = 1, npaire = 1
     aam211t=aamt; aam211t(ny-2:ny,:)=aam11t(ny-2:ny,:)
     bbm211t=bbmt; bbm211t(ny-2:ny,:)=bbm11t(ny-2:ny,:)
     ccm211t=ccmt; ccm211t(ny-2:ny,:)=ccm11t(ny-2:ny,:)
     ddm211t=ddmt; ddm211t(ny-2:ny,:)=ddm11t(ny-2:ny,:)
     eem211t=eemt; eem211t(ny-2:ny,:)=eem11t(ny-2:ny,:)
     qqm211t=qqmt; qqm211t(ny-2:ny,:)=qqm11t(ny-2:ny,:)
     rrm211t=rrmt; rrm211t(ny-2:ny,:)=rrm11t(ny-2:ny,:)

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
     aam(ny-2  )=-two*(asty+bsty)
     aam(4     )=-two*(as4y+bs4y+cs4y)
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
     ccm(2     )=zero!bs2y
     ccm(ny-1  )=zero!bsmy
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
     ! velocity, ncly1 = 2, nclyn = 2
     call ludecomp7(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,&
          vvm,wwm,zzm,ny)
     ! velocity, ncly1 = 1, nclyn = 1, npaire = 0
     call ludecomp7(aam10,bbm10,ccm10,ddm10,eem10,qqm10,ggm10,hhm10,ssm10,rrm10,&
          vvm10,wwm10,zzm10,ny)
     ! velocity, ncly1 = 1, nclyn = 1, npaire = 1
     call ludecomp7(aam11,bbm11,ccm11,ddm11,eem11,qqm11,ggm11,hhm11,ssm11,rrm11,&
          vvm11,wwm11,zzm11,ny)
     ! velocity, ncly1 = 0, nclyn = 0
     call ludecomp7(aam0,bbm0,ccm0,ddm0,eem0,qqm0,ggm0,hhm0,ssm0,rrm0,&
          vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ny)
     ! velocity, ncly1 = 1, nclyn = 2, npaire = 0
     call ludecomp7(aam120,bbm120,ccm120,ddm120,eem120,qqm120,ggm120,hhm120,ssm120,rrm120,&
          vvm120,wwm120,zzm120,ny)
     ! velocity, ncly1 = 1, nclyn = 2, npaire = 1
     call ludecomp7(aam121,bbm121,ccm121,ddm121,eem121,qqm121,ggm121,hhm121,ssm121,rrm121,&
          vvm121,wwm121,zzm121,ny)
     ! velocity, ncly1 = 2, nclyn = 1, npaire = 0
     call ludecomp7(aam210,bbm210,ccm210,ddm210,eem210,qqm210,ggm210,hhm210,ssm210,rrm210,&
          vvm210,wwm210,zzm210,ny)
     ! velocity, ncly1 = 2, nclyn = 1, npaire = 1
     call ludecomp7(aam211,bbm211,ccm211,ddm211,eem211,qqm211,ggm211,hhm211,ssm211,rrm211,&
          vvm211,wwm211,zzm211,ny)
     ! Scalars
     do is = 1, numscalar
        ! scalar, ncly1 = 2, nclyn = 2
        call ludecomp7(aamt(:,is),bbmt(:,is),ccmt(:,is),ddmt(:,is),eemt(:,is),&
             qqmt(:,is),ggmt(:,is),hhmt(:,is),ssmt(:,is),rrmt(:,is),&
             vvmt(:,is),wwmt(:,is),zzmt(:,is),ny)
        ! scalar, ncly1 = 1, nclyn = 1, npaire = 0
        call ludecomp7(aam10t(:,is),bbm10t(:,is),ccm10t(:,is),ddm10t(:,is),eem10t(:,is),&
             qqm10t(:,is),ggm10t(:,is),hhm10t(:,is),ssm10t(:,is),rrm10t(:,is),&
             vvm10t(:,is),wwm10t(:,is),zzm10t(:,is),ny)
        ! scalar, ncly1 = 1, nclyn = 1, npaire = 1
        call ludecomp7(aam11t(:,is),bbm11t(:,is),ccm11t(:,is),ddm11t(:,is),eem11t(:,is),&
             qqm11t(:,is),ggm11t(:,is),hhm11t(:,is),ssm11t(:,is),rrm11t(:,is),&
             vvm11t(:,is),wwm11t(:,is),zzm11t(:,is),ny)
        ! scalar, ncly1 = 0, nclyn = 0
        call ludecomp7(aam0t(:,is),bbm0t(:,is),ccm0t(:,is),ddm0t(:,is),eem0t(:,is),&
             qqm0t(:,is),ggm0t(:,is),hhm0t(:,is),ssm0t(:,is),rrm0t(:,is),&
             vvm0t(:,is),wwm0t(:,is),zzm0t(:,is),l1mt(:,is),l2mt(:,is),l3mt(:,is),u1mt(:,is),u2mt(:,is),u3mt(:,is),ny)
        ! scalar, ncly1 = 1, nclyn = 2, npaire = 0
        call ludecomp7(aam120t(:,is),bbm120t(:,is),ccm120t(:,is),ddm120t(:,is),eem120t(:,is),&
             qqm120t(:,is),ggm120t(:,is),hhm120t(:,is),ssm120t(:,is),rrm120t(:,is),&
             vvm120t(:,is),wwm120t(:,is),zzm120t(:,is),ny)
        ! scalar, ncly1 = 1, nclyn = 2, npaire = 1
        call ludecomp7(aam121t(:,is),bbm121t(:,is),ccm121t(:,is),ddm121t(:,is),eem121t(:,is),&
             qqm121t(:,is),ggm121t(:,is),hhm121t(:,is),ssm121t(:,is),rrm121t(:,is),&
             vvm121t(:,is),wwm121t(:,is),zzm121t(:,is),ny)
        ! scalar, ncly1 = 2, nclyn = 1, npaire = 0
        call ludecomp7(aam210t(:,is),bbm210t(:,is),ccm210t(:,is),ddm210t(:,is),eem210t(:,is),&
             qqm210t(:,is),ggm210t(:,is),hhm210t(:,is),ssm210t(:,is),rrm210t(:,is),&
             vvm210t(:,is),wwm210t(:,is),zzm210t(:,is),ny)
        ! scalar, ncly1 = 2, nclyn = 1, npaire = 1
        call ludecomp7(aam211t(:,is),bbm211t(:,is),ccm211t(:,is),ddm211t(:,is),eem211t(:,is),&
             qqm211t(:,is),ggm211t(:,is),hhm211t(:,is),ssm211t(:,is),rrm211t(:,is),&
             vvm211t(:,is),wwm211t(:,is),zzm211t(:,is),ny)
     enddo
  else
     call ludecomp9(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ttm,uum,sssm,zzzm,ny)
     !NEED TO BE DONE: deal with other cases
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! END MATRIX M2 TIME IMPLICIT    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine implicit_schemes

!
! Allocate 1D arrays containing LU decompositions
!
subroutine init_implicit()

  USE decomp_2d
  USE param
  USE variables
  implicit none

  ! velocity, ncly1 = 2, nclyn = 2
  allocate(aam(ny),bbm(ny),ccm(ny),ddm(ny),eem(ny),ggm(ny),hhm(ny),wwm(ny),zzm(ny))
  allocate(rrm(ny),qqm(ny),vvm(ny),ssm(ny))
  allocate(sssm(ny),zzzm(ny),ttm(ny),uum(ny)) ! nona
  ! velocity, ncly1 = 1, nclyn = 1, npaire = 0
  allocate(aam10(ny),bbm10(ny),ccm10(ny),ddm10(ny),eem10(ny),ggm10(ny),hhm10(ny),wwm10(ny),zzm10(ny))
  allocate(rrm10(ny),qqm10(ny),vvm10(ny),ssm10(ny))
  ! velocity, ncly1 = 1, nclyn = 1, npaire = 1
  allocate(aam11(ny),bbm11(ny),ccm11(ny),ddm11(ny),eem11(ny),ggm11(ny),hhm11(ny),wwm11(ny),zzm11(ny))
  allocate(rrm11(ny),qqm11(ny),vvm11(ny),ssm11(ny))
  ! velocity, ncly1 = 0, nclyn = 0
  allocate(aam0(ny),bbm0(ny),ccm0(ny),ddm0(ny),eem0(ny),ggm0(ny),hhm0(ny),wwm0(ny),zzm0(ny))
  allocate(rrm0(ny),qqm0(ny),vvm0(ny),ssm0(ny),l1m(ny),l2m(ny),l3m(ny),u1m(ny),u2m(ny),u3m(ny))
  ! velocity, ncly1 = 1, nclyn = 2, npaire = 0
  allocate(aam120(ny),bbm120(ny),ccm120(ny),ddm120(ny),eem120(ny),ggm120(ny),hhm120(ny),wwm120(ny),zzm120(ny))
  allocate(rrm120(ny),qqm120(ny),vvm120(ny),ssm120(ny))
  ! velocity, ncly1 = 1, nclyn = 2, npaire = 1
  allocate(aam121(ny),bbm121(ny),ccm121(ny),ddm121(ny),eem121(ny),ggm121(ny),hhm121(ny),wwm121(ny),zzm121(ny))
  allocate(rrm121(ny),qqm121(ny),vvm121(ny),ssm121(ny))
  ! velocity, ncly1 = 2, nclyn = 1, npaire = 0
  allocate(aam210(ny),bbm210(ny),ccm210(ny),ddm210(ny),eem210(ny),ggm210(ny),hhm210(ny),wwm210(ny),zzm210(ny))
  allocate(rrm210(ny),qqm210(ny),vvm210(ny),ssm210(ny))
  ! velocity, ncly1 = 2, nclyn = 1, npaire = 1
  allocate(aam211(ny),bbm211(ny),ccm211(ny),ddm211(ny),eem211(ny),ggm211(ny),hhm211(ny),wwm211(ny),zzm211(ny))
  allocate(rrm211(ny),qqm211(ny),vvm211(ny),ssm211(ny))
  ! scalar, ncly1 = 2, nclyn = 2
  allocate(aamt(ny,numscalar),bbmt(ny,numscalar),ccmt(ny,numscalar),ddmt(ny,numscalar),eemt(ny,numscalar))
  allocate(ggmt(ny,numscalar),hhmt(ny,numscalar),wwmt(ny,numscalar),zzmt(ny,numscalar))
  allocate(rrmt(ny,numscalar),qqmt(ny,numscalar),vvmt(ny,numscalar),ssmt(ny,numscalar))
  allocate(uumt(ny,numscalar),ttmt(ny,numscalar),sssmt(ny,numscalar),zzzmt(ny,numscalar)) ! nona
  ! scalar, ncly1 = 1, nclyn = 1, npaire = 0
  allocate(aam10t(ny,numscalar),bbm10t(ny,numscalar),ccm10t(ny,numscalar),ddm10t(ny,numscalar),eem10t(ny,numscalar))
  allocate(ggm10t(ny,numscalar),hhm10t(ny,numscalar),wwm10t(ny,numscalar),zzm10t(ny,numscalar))
  allocate(rrm10t(ny,numscalar),qqm10t(ny,numscalar),vvm10t(ny,numscalar),ssm10t(ny,numscalar))
  ! scalar, ncly1 = 1, nclyn = 1, npaire = 1
  allocate(aam11t(ny,numscalar),bbm11t(ny,numscalar),ccm11t(ny,numscalar),ddm11t(ny,numscalar),eem11t(ny,numscalar))
  allocate(ggm11t(ny,numscalar),hhm11t(ny,numscalar),wwm11t(ny,numscalar),zzm11t(ny,numscalar))
  allocate(rrm11t(ny,numscalar),qqm11t(ny,numscalar),vvm11t(ny,numscalar),ssm11t(ny,numscalar))
  ! scalar, ncly1 = 0, nclyn = 0
  allocate(aam0t(ny,numscalar),bbm0t(ny,numscalar),ccm0t(ny,numscalar),ddm0t(ny,numscalar),eem0t(ny,numscalar))
  allocate(ggm0t(ny,numscalar),hhm0t(ny,numscalar),wwm0t(ny,numscalar),zzm0t(ny,numscalar))
  allocate(rrm0t(ny,numscalar),qqm0t(ny,numscalar),vvm0t(ny,numscalar),ssm0t(ny,numscalar))
  allocate(l1mt(ny,numscalar),l2mt(ny,numscalar),l3mt(ny,numscalar),u1mt(ny,numscalar),u2mt(ny,numscalar),u3mt(ny,numscalar))
  ! scalar, ncly1 = 1, nclyn = 2, npaire = 0
  allocate(aam120t(ny,numscalar),bbm120t(ny,numscalar),ccm120t(ny,numscalar),ddm120t(ny,numscalar),eem120t(ny,numscalar))
  allocate(ggm120t(ny,numscalar),hhm120t(ny,numscalar),wwm120t(ny,numscalar),zzm120t(ny,numscalar))
  allocate(rrm120t(ny,numscalar),qqm120t(ny,numscalar),vvm120t(ny,numscalar),ssm120t(ny,numscalar))
  ! scalar, ncly1 = 1, nclyn = 2, npaire = 1
  allocate(aam121t(ny,numscalar),bbm121t(ny,numscalar),ccm121t(ny,numscalar),ddm121t(ny,numscalar),eem121t(ny,numscalar))
  allocate(ggm121t(ny,numscalar),hhm121t(ny,numscalar),wwm121t(ny,numscalar),zzm121t(ny,numscalar))
  allocate(rrm121t(ny,numscalar),qqm121t(ny,numscalar),vvm121t(ny,numscalar),ssm121t(ny,numscalar))
  ! scalar, ncly1 = 2, nclyn = 1, npaire = 0
  allocate(aam210t(ny,numscalar),bbm210t(ny,numscalar),ccm210t(ny,numscalar),ddm210t(ny,numscalar),eem210t(ny,numscalar))
  allocate(ggm210t(ny,numscalar),hhm210t(ny,numscalar),wwm210t(ny,numscalar),zzm210t(ny,numscalar))
  allocate(rrm210t(ny,numscalar),qqm210t(ny,numscalar),vvm210t(ny,numscalar),ssm210t(ny,numscalar))
  ! scalar, ncly1 = 2, nclyn = 1, npaire = 1
  allocate(aam211t(ny,numscalar),bbm211t(ny,numscalar),ccm211t(ny,numscalar),ddm211t(ny,numscalar),eem211t(ny,numscalar))
  allocate(ggm211t(ny,numscalar),hhm211t(ny,numscalar),wwm211t(ny,numscalar),zzm211t(ny,numscalar))
  allocate(rrm211t(ny,numscalar),qqm211t(ny,numscalar),vvm211t(ny,numscalar),ssm211t(ny,numscalar))

end subroutine init_implicit

!
! Used to build the scalar implicit coefficients
!
subroutine init_implicit_coef(tab1d, tab2d)

  use decomp_2d, only : mytype
  use variables, only : ny, numscalar

  implicit none

  real(mytype), dimension(ny), intent(in) :: tab1d
  real(mytype), dimension(ny,numscalar), intent(out) :: tab2d

  integer :: is

  do is = 1, numscalar
     tab2d(:,is) = tab1d(:)
  enddo

end subroutine init_implicit_coef

end module ydiff_implicit
