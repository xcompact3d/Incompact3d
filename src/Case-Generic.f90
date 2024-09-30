!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module case_generic

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  real(mytype), save, allocatable, dimension(:,:,:) :: vol1,volSimps1
  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_generic, &
            boundary_conditions_generic, &
            postprocess_generic, &
            visu_generic

contains

  subroutine init_generic (ux1,uy1,uz1,ep1,phi1)

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
  end subroutine init_generic

  subroutine boundary_conditions_generic (ux,uy,uz,phi,ep)

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

  end subroutine boundary_conditions_generic

  subroutine postprocess_generic(ux1,uy1,uz1,phi1,ep1)

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

  end subroutine postprocess_generic

  !############################################################################
  !!
  !!  SUBROUTINE: visu_generic
  !!      AUTHOR: CF
  !! DESCRIPTION: Performs case-specific visualization
  !!
  !############################################################################
  subroutine visu_generic(ux1, uy1, uz1, pp3, phi1, ep1, num)

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
    integer, intent(in) :: num

  end subroutine visu_generic

end module case_generic
