!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module pipe

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  use variables
  use param
  !!$use complex_geometry, only: tol

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: geomcomplex_pipe, init_pipe, boundary_conditions_pipe, postprocess_pipe

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: geomcomplex_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Generates epsilon matrix for pipe geometry
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine geomcomplex_pipe(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,dz,remp)
  !
  !********************************************************************

    use MPI
    use param,only : zero,one,two,yly,zlz
    use ibm_param

    implicit none

    integer                                         :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny)                      :: yp
    real(mytype)                                    :: dx,dz,remp,tol
    !LOCALS
    real(mytype)                                    :: r,ym,zm,yc,zc
    integer                                         :: i,j,k,code,ierror

    epsi(:,:,:) = zero
    yc = yly/two
    zc = zlz/two
    tol=1e-15

    !safety check
    if (nrank == 0) then
       if (rai.le.0) then
           write(*,*) 'SIMULATION IS STOPPED!'
           write(*,*) 'Please specify a valid value for the pipe inner radius in input.i3d (rai)'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif
       if (rao.le.0) then
           write(*,*) 'SIMULATION IS STOPPED!'
           write(*,*) 'Please specify a valid value for the pipe outer radius in input.i3d (rao)'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif
    endif

    do k=nzi,nzf
        zm=real(k-1,mytype)*dz-zc
        do j=nyi,nyf
            ym=yp(j)-yc
            r=sqrt(ym*ym+zm*zm)
            do i=nxi,nxf
                if (r.gt.rai.and.r.lt.rao) then
                   epsi(i,j,k)=remp
                elseif (abs(r-rai).lt.tol) then
                    epsi(i,j,k)=remp
                elseif (abs(r-rao).lt.tol) then
                    epsi(i,j,k)=remp
                endif
            enddo
        enddo
    enddo
    !
    return

  end subroutine geomcomplex_pipe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  subroutine: init_pipe
  !!      AUTHOR: Rodrigo Vicente Cruz
  !! DESCRIPTION: Initial laminar conditions in pipe
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !********************************************************************
  !
  subroutine init_pipe (ux1,uy1,uz1,ep1,phi1)
  !
  !********************************************************************

    use decomp_2d_io
    use variables
    use param
    use ibm_param
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3))              :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar)    :: phi1

    real(mytype)                    :: r,ym,zm,theta
    real(mytype)                    :: um,yc,zc
    real(mytype)                    :: nu
    integer                         :: k,j,i,fh,ierror,ii,is,it,icht,code
    integer (kind=MPI_OFFSET_KIND)  :: disp
    !
    yc = yly / two
    zc = zlz / two

    if (iscalar.ne.0) then
        !Analytical laminar temperature profile, with nu=4.36
        nu=real(48./11.,8)
        do is=1,numscalar
            do k=1,xsize(3) !MIXED-TYPE BOUNDARY CONDITION : PHI = (Tw-T)/(Tw-Tb)
                zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
                do j=1,xsize(2)
                    if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
                    if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
                    r=sqrt(ym*ym+zm*zm)
                    do i=1,xsize(1)
                        if (r.le.rai.and.ep1(i,j,k).eq.0) then
                            phi1(i,j,k,is)=two*nu*(three/sixteen + r**four - r**two)
                        else
                            phi1(i,j,k,is)=zero
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

    ux1=zero;uy1=zero;uz1=zero
    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))
       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
    endif
    !modulation of the random noise + initial velocity profile
    do k=1,xsize(3)
        zm=dz*real(xstart(3)-1+k-1,mytype)-zc 
        do j=1,xsize(2)
            if (istret.eq.0) ym=real(j+xstart(2)-1-1,mytype)*dy-yc
            if (istret.ne.0) ym=yp(j+xstart(2)-1)-yc
            r=sqrt(ym*ym+zm*zm)    
            um=exp(-fifteen*r*r)
            !Poiseuille flow
            bxx1(j,k)=two*(one-(ym**two+zm**two)/(rai**two))
            bxy1(j,k)=zero
            bxz1(j,k)=zero
            do i=1,xsize(1)
                if (r.le.rai.and.ep1(i,j,k).eq.0) then
                    ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+bxx1(j,k)
                    uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)+bxy1(j,k)
                    uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)+bxz1(j,k)
                else
                    ux1(i,j,k)=zero
                    uy1(i,j,k)=zero
                    uz1(i,j,k)=zero
                endif
            enddo
        enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif
    return

  end subroutine init_pipe
  !********************************************************************
  !
  subroutine boundary_conditions_pipe (ux,uy,uz,phi)
  !
  !********************************************************************

    use param
    use variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    !
    return
    !
  end subroutine boundary_conditions_pipe
  !********************************************************************
  !
  subroutine postprocess_pipe(ux1,uy1,uz1,pp3,phi1,ep1)
  !
  !********************************************************************

    use var, ONLY : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_pipe

end module pipe
