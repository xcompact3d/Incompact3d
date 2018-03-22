!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contain duplicated code that scatters data from the
! MPI_ALLTOALLV receive buffer to destinations. It is 'included' by two 
! subroutines in decomp_2d.f90

! Note:
!   in  --> receive buffer
!   out --> destination array
!   pos --> pointer for the receive buffer
!    - for normal ALLTOALLV, points to the beginning of receive buffer (=1)
!    - for shared memory code, note the receive buffer is shared by all cores
!      on same node, so 'pos' points to the correct location for this core

    integer, intent(IN) :: ndir
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos
    
#ifndef SHM
    pos = 1
#endif
    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       endif

       if (ndir==1) then
#ifdef SHM
          pos = decomp%y1disp_o(m) + 1
#endif
          do k=1,n3
             do j=i1,i2
                do i=1,n1
                   out(i,j,k) = in(pos)
                   pos = pos + 1
                enddo
             enddo
          enddo
       else if (ndir==2) then
#ifdef SHM
          pos = decomp%z2disp_o(m) + 1
#endif
          do k=i1,i2
             do j=1,n2
                do i=1,n1
                   out(i,j,k) = in(pos)
                   pos = pos + 1
                enddo
             enddo
          enddo
       else if (ndir==3) then
#ifdef SHM
          pos = decomp%y2disp_o(m) + 1
#endif
          do k=1,n3
             do j=i1,i2
                do i=1,n1
                   out(i,j,k) = in(pos)
                   pos = pos + 1
                enddo
             enddo
          enddo
       else if (ndir==4) then
#ifdef SHM
          pos = decomp%x1disp_o(m) + 1
#endif
          do k=1,n3
             do j=1,n2
                do i=i1,i2
                   out(i,j,k) = in(pos)
                   pos = pos + 1
                enddo
             enddo
          enddo
       endif
    enddo
