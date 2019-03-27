!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!        FILE: BC-Jet.f90
!!!      AUTHOR: Arash
!!!    MODIFIED: Paul Bartholomew
!!! DESCRIPTION: This module describes the jet case.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module jet

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  !probes
  integer, save :: nprobes, ntimes1, ntimes2
  integer, save, allocatable, dimension(:) :: rankprobes, nxprobes, nyprobes, nzprobes

  real(mytype),save,allocatable,dimension(:) :: usum,vsum,wsum,uusum,uvsum,uwsum,vvsum,vwsum,wwsum

  REAL(mytype) :: outflow

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_jet, boundary_conditions_jet, postprocessing_jet, momentum_forcing_jet

contains

  subroutine init_jet (rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,phis1,phiss1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE var, ONLY : ta2, tb2, tc2, td2

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1,phis1,phiss1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1,drho1

    real(mytype) :: y
    integer :: i, j, k, ii, is, code
    real(mytype) :: um

    integer, dimension (:), allocatable :: seed

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
       do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)
          um=exp(-y**2)
          um = one
          do i=1,xsize(1)
             ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)
             uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
          enddo
       enddo
    enddo

    call boundary_conditions_jet (rho1,ux1,uy1,uz1,phi1)

    !INIT FOR G AND U=MEAN FLOW + NOISE
    if (xstart(2)==1) then
      j = 1
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          ux1(i, j, k) = byx1(i, k)
          uy1(i, j, k) = byy1(i, k)
          uz1(i, j, k) = byz1(i, k)
        enddo
      enddo
    endif
    call transpose_x_to_y(ux1, ta2)
    call transpose_x_to_y(uy1, tb2)
    call transpose_x_to_y(uz1, tc2)
    call transpose_x_to_y(rho1(:,:,:,1), td2)
    do k=1,ysize(3)
       do j=2,ysize(2)
          do i=1,ysize(1)
             ta2(i,j,k)=ta2(i,j,k)+ta2(i,1,k)
             tb2(i,j,k)=tb2(i,j,k)+tb2(i,1,k)
             tc2(i,j,k)=tc2(i,j,k)+tc2(i,1,k)

             td2(i,j,k)=td2(i,1,k)
          enddo
       enddo
    enddo
    call transpose_y_to_x(ta2, ux1)
    call transpose_y_to_x(tb2, uy1)
    call transpose_y_to_x(tc2, uz1)
    call transpose_y_to_x(td2, rho1(:,:,:,1))
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             dux1(i,j,k,1)=ux1(i,j,k)
             duy1(i,j,k,1)=uy1(i,j,k)
             duz1(i,j,k,1)=uz1(i,j,k)

             drho1(i,j,k,1) = rho1(i,j,k,1)
             do is = 2, ntime
                dux1(i,j,k,is)=dux1(i,j,k,is - 1)
                duy1(i,j,k,is)=duy1(i,j,k,is - 1)
                duz1(i,j,k,is)=duz1(i,j,k,is - 1)

                drho1(i,j,k,is) = drho1(i,j,k,is - 1)
             enddo
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_jet
  !********************************************************************
  subroutine boundary_conditions_jet (rho,ux,uy,uz,phi)

    USE MPI
    USE param
    USE variables
    USE decomp_2d
    !USE var

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho


    integer :: i, j, k
    real(mytype) :: D, r, x, z, r1, r2, x1, z1
    real(mytype) :: inflow, perturbation
    real(mytype) :: timeswitch

    integer :: ierr
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods

    real(mytype) :: uu1,uv1,uw1,x2,y1,y2,ya,y,xc,zc,yc

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, ierr)


    D = one
    perturbation = zero

    if (t.lt.one) then
       timeswitch = sin(t * pi / two)
    else
       timeswitch = one
    endif
    timeswitch = one !! Nichols starts with a jet 'column'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Set inflow
    inflow = zero
    j = 1
    if (xstart(2) == 1) then
      do k = 1, xsize(3)
        z = real(k + xstart(3) - 2, mytype) * dz - half * zlz
        do i = 1, xsize(1)
          x = real(i + xstart(1) - 2, mytype) * dx - half * xlx
          r = sqrt(x**2 + z**2) + 1.e-16_mytype !! Can't divide by zero!

          !! Set mean inflow
          byx1(i, k) = zero
          byy1(i, k) = (u1 - u2) * half &
               * (one + tanh((12.5_mytype / four) * ((D / two) / r - two * r / D))) + u2
          byz1(i, k) = zero

          rho(i, 1, k, 1) = (dens1 - dens2) * half &
               * (one + tanh((12.5_mytype / four) * ((D / two) / r - two * r / D))) + dens2

          !! Apply transient behaviour
          ! if (r.lt.D/two) then
          !   perturbation = inflow_noise * sin(r * x * z * t)
          ! else
          !   perturbation = zero
          ! endif
          ! byy1(i, k) = byy1(i, k) + perturbation
          byy1(i, k) = timeswitch * byy1(i, k)

          !! Sum up total flux
          inflow = inflow + byy1(i, k)
        enddo
      enddo
    endif

    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
    bxo(:,:) = two * bxo(:,:) - one
    byo(:,:) = two * byo(:,:) - one
    bzo(:,:) = two * bzo(:,:) - one
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxx1(j,k)=bxo(j,k)*inflow_noise
        bxy1(j,k)=byo(j,k)*inflow_noise
        bxz1(j,k)=bzo(j,k)*inflow_noise
      enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Apply lateral boundary conditions
    !! XXX Assume purely radial flow
    !! XXX In X-pencils
    
    xc=half*xlx
    zc=half*xlx

    !! X-BC
    IF (nclx1.EQ.2) THEN
       if(xstart(1).eq.1)then!
          x = -xc
          i = 1
          do k=1,xsize(3)
             z=real(k + xstart(3) - 2,mytype)*dz-zc
             x2=x
             y2=z
             x1=x+dx
             y1=y2*x1/x2
             r1=sqrt(x1**2+y1**2)
             r2=sqrt(x2**2+y2**2)
             if(r1.gt.r2)print*,'bug CL'
             if(k.eq.1)then!cas premier point
                do j=1,xsize(2)

                   bxx1(j,k)=r1*ux(i + 1,j,k+1)/r2
                   bxy1(j,k)=   uy(i + 1,j,k+1)
                   bxz1(j,k)=r1*uz(i + 1,j,k+1)/r2
                enddo
             elseif(k.eq.(nz-1)/2+1)then!cas point du milieu
                do j=1,xsize(2)

                   bxx1(j,k)=r1*ux(i+1,j,k)/r2
                   bxy1(j,k)=   uy(i+1,j,k)
                   bxz1(j,k)=r1*uz(i+1,j,k)/r2
                enddo
             elseif(k.eq.nx)then!cas dernier point
                do j=1,xsize(2)

                   bxx1(j,k)=r1*ux(i+1,j,k-1)/r2
                   bxy1(j,k)=   uy(i+1,j,k-1)
                   bxz1(j,k)=r1*uz(i+1,j,k-1)/r2
                enddo
             else!cas general
                if    (z.gt.0.)then
                   ya=y2-dz
                   do j=1,xsize(2)
                      uu1=(ux(i+1,j,k)-ux(i+1,j,k-1))*(y1-ya)/(y2-ya)+ux(i+1,j,k-1)
                      uv1=(uy(i+1,j,k)-uy(i+1,j,k-1))*(y1-ya)/(y2-ya)+uy(i+1,j,k-1)
                      uw1=(uz(i+1,j,k)-uz(i+1,j,k-1))*(y1-ya)/(y2-ya)+uz(i+1,j,k-1)

                      bxx1(j,k)=r1*uu1/r2
                      bxy1(j,k)=   uv1
                      bxz1(j,k)=r1*uw1/r2
                   enddo
                elseif(z.lt.0.)then
                   ya=y2+dz
                   do j=1,xsize(2)
                      uu1=(ux(i+1,j,k+1)-ux(2,j,k))*(y1-ya)/(ya-y2)+ux(i+1,j,k+1)
                      uv1=(uy(i+1,j,k+1)-uy(2,j,k))*(y1-ya)/(ya-y2)+uy(i+1,j,k+1)
                      uw1=(uz(i+1,j,k+1)-uz(2,j,k))*(y1-ya)/(ya-y2)+uz(i+1,j,k+1)

                      bxx1(j,k)=r1*uu1/r2
                      bxy1(j,k)=   uv1
                      bxz1(j,k)=r1*uw1/r2
                   enddo
                endif
             endif
          enddo
       endif
    ENDIF

    IF (nclxn.EQ.2) THEN
       if(xend(1).eq.nx)then
          x=xc
          i = xsize(1)
          do k=1,xsize(3)
             z=real(k + xstart(3) - 2)*dz-zc
             x2=x
             y2=z
             x1=x-dx
             y1=y2*x1/x2
             r1=sqrt(x1**2+y1**2)
             r2=sqrt(x2**2+y2**2)
             if(r1.gt.r2)print*,'bug CL'
             if(k.eq.1)then!cas premier point
                do j=1,xsize(2)

                   bxxn(j,k)=r1*ux(i-1,j,k+1)/r2
                   bxyn(j,k)=   uy(i-1,j,k+1)
                   bxzn(j,k)=r1*uz(i-1,j,k+1)/r2
                enddo
             elseif(k.eq.(nz-1)/2+1)then!cas point du milieu
                do j=1,xsize(2)

                   bxxn(j,k)=r1*ux(i-1,j,k)/r2
                   bxyn(j,k)=   uy(i-1,j,k)
                   bxzn(j,k)=r1*uz(i-1,j,k)/r2
                enddo
             elseif(k.eq.nz)then!cas dernier point
                do j=1,xsize(2)

                   bxxn(j,k)=r1*ux(i-1,j,k-1)/r2
                   bxyn(j,k)=   uy(i-1,j,k-1)
                   bxzn(j,k)=r1*uz(i-1,j,k-1)/r2
                enddo
             else !cas general
                if (z.gt.0.) then
                   ya=y2-dz
                   do j=1,xsize(2)
                      uu1=(ux(i-1,j,k)-ux(i-1,j,k-1))*(y1-ya)/(y2-ya)+ux(i-1,j,k-1)
                      uv1=(uy(i-1,j,k)-uy(i-1,j,k-1))*(y1-ya)/(y2-ya)+uy(i-1,j,k-1)
                      uw1=(uz(i-1,j,k)-uz(i-1,j,k-1))*(y1-ya)/(y2-ya)+uz(i-1,j,k-1)

                      bxxn(j,k)=r1*uu1/r2
                      bxyn(j,k)=   uv1
                      bxzn(j,k)=r1*uw1/r2
                   enddo
                elseif(z.lt.0.)then
                   ya=y2+dz
                   do j=1,xsize(2)
                      uu1=(ux(i-1,j,k+1)-ux(i-1,j,k))*(y1-ya)/(ya-y2)+ux(i-1,j,k+1)
                      uv1=(uy(i-1,j,k+1)-uy(i-1,j,k))*(y1-ya)/(ya-y2)+uy(i-1,j,k+1)
                      uw1=(uz(i-1,j,k+1)-uz(i-1,j,k))*(y1-ya)/(ya-y2)+uz(i-1,j,k+1)

                      bxxn(j,k)=r1*uu1/r2
                      bxyn(j,k)=   uv1
                      bxzn(j,k)=r1*uw1/r2
                   enddo
                endif
             endif
          enddo
       endif
    ENDIF

    !! Z-BC
    IF ((nclz1.EQ.2).AND.(xstart(3).EQ.1)) THEN
       k = 1
       z = -zc
       do i=1,xsize(1)
          x=real(i + xstart(1) - 2, mytype) * dx - xc
          x2=z
          y2=x
          x1=z+dz
          y1=y2*x1/x2
          r1=sqrt(x1**2+y1**2)
          r2=sqrt(x2**2+y2**2)
          if(r1.gt.r2)print*,'bug CL'
          if(i.eq.1)then!cas premier point
             do j=1,xsize(2)

                bzx1(i,j)=r1*ux(i+1,j,k + 1)/r2
                bzy1(i,j)=   uy(i+1,j,k + 1)
                bzz1(i,j)=r1*uz(i+1,j,k + 1)/r2
             enddo
          elseif(i.eq.(nx-1)/2+1)then!cas point du milieu
             do j=1,xsize(2)

                bzx1(i,j)=r1*ux(i,j,k + 1)/r2
                bzy1(i,j)=   uy(i,j,k + 1)
                bzz1(i,j)=r1*uz(i,j,k + 1)/r2
             enddo
          elseif(i.eq.nx)then!cas dernier point
             do j=1,xsize(2)

                bzx1(i,j)=r1*ux(i-1,j,k + 1)/r2
                bzy1(i,j)=   uy(i-1,j,k + 1)
                bzz1(i,j)=r1*uz(i-1,j,k + 1)/r2
             enddo
          else!cas general
             if    (x.gt.0.)then
                ya=y2-dx
                do j=1,xsize(2)

                   uu1=(ux(i,j,k + 1)-ux(i-1,j,k + 1))*(y1-ya)/(y2-ya)+ux(i-1,j,k + 1)
                   uv1=(uy(i,j,k + 1)-uy(i-1,j,k + 1))*(y1-ya)/(y2-ya)+uy(i-1,j,k + 1)
                   uw1=(uz(i,j,k + 1)-uz(i-1,j,k + 1))*(y1-ya)/(y2-ya)+uz(i-1,j,k + 1)

                   bzx1(i,j)=r1*uu1/r2
                   bzy1(i,j)=   uv1
                   bzz1(i,j)=r1*uw1/r2
                enddo
             elseif(x.lt.0.)then
                ya=y2+dx
                do j=1,xsize(2)

                   uu1=(ux(i+1,j,k + 1)-ux(i,j,k + 1))*(y1-ya)/(ya-y2)+ux(i+1,j,k + 1)
                   uv1=(uy(i+1,j,k + 1)-uy(i,j,k + 1))*(y1-ya)/(ya-y2)+uy(i+1,j,k + 1)
                   uw1=(uz(i+1,j,k + 1)-uz(i,j,k + 1))*(y1-ya)/(ya-y2)+uz(i+1,j,k + 1)

                   bzx1(i,j)=r1*uu1/r2
                   bzy1(i,j)=   uv1
                   bzz1(i,j)=r1*uw1/r2
                enddo
             endif
          endif
       enddo
    ENDIF

    IF ((nclzn.EQ.2).AND.(xend(3).EQ.nz)) THEN
       z=zc
       k = xsize(3)
       do i=1,xsize(1)
          x=real(i + xstart(1) - 2, mytype)*dx-xc
          x2=z
          y2=x
          x1=z-dz
          y1=y2*x1/x2
          r1=sqrt(x1**2+y1**2)
          r2=sqrt(x2**2+y2**2)
          if(r1.gt.r2)print*,'bug CL'
          if(i.eq.1)then!cas premier point
             do j=1,xsize(2)

                bzxn(i,j)=r1*ux(i+1,j,k-1)/r2
                bzyn(i,j)=   uy(i+1,j,k-1)
                bxzn(i,j)=r1*uz(i+1,j,k-1)/r2
             enddo
          elseif(i.eq.(nx-1)/2+1)then!cas point du milieu
             do j=1,xsize(2)

                bzxn(i,j)=r1*ux(i,j,k-1)/r2
                bzyn(i,j)=   uy(i,j,k-1)
                bzzn(i,j)=r1*uz(i,j,k-1)/r2
             enddo
          elseif(i.eq.nx)then!cas dernier point
             do j=1,xsize(2)

                bzxn(i,k)=r1*ux(i-1,j,k-1)/r2
                bzyn(i,k)=   uy(i-1,j,k-1)
                bzzn(i,k)=r1*uz(i-1,j,k-1)/r2
             enddo
          else !cas general
             if (x.gt.0.) then
                ya=y2-dx
                do j=1,xsize(2)

                   uu1=(ux(i,j,k-1)-ux(i-1,j,k-1))*(y1-ya)/(y2-ya)+ux(i-1,j,k-1)
                   uv1=(uy(i,j,k-1)-uy(i-1,j,k-1))*(y1-ya)/(y2-ya)+uy(i-1,j,k-1)
                   uw1=(uz(i,j,k-1)-uz(i-1,j,k-1))*(y1-ya)/(y2-ya)+uz(i-1,j,k-1)

                   bzxn(i,k)=r1*uu1/r2
                   bzyn(i,k)=   uv1
                   bzzn(i,k)=r1*uw1/r2
                enddo
             elseif(x.lt.0.)then
                ya=y2+dx
                do j=1,xsize(2)
                   uu1=(ux(i-1,j,k+1)-ux(i-1,j,k))*(y1-ya)/(ya-y2)+ux(i-1,j,k+1)
                   uv1=(uy(i-1,j,k+1)-uy(i-1,j,k))*(y1-ya)/(ya-y2)+uy(i-1,j,k+1)
                   uw1=(uz(i-1,j,k+1)-uz(i-1,j,k))*(y1-ya)/(ya-y2)+uz(i-1,j,k+1)

                   uu1=(ux(i+1,j,k-1)-ux(i,j,k-1))*(y1-ya)/(ya-y2)+ux(i+1,j,k-1)
                   uv1=(uy(i+1,j,k-1)-uy(i,j,k-1))*(y1-ya)/(ya-y2)+uy(i+1,j,k-1)
                   uw1=(uz(i+1,j,k-1)-uz(i,j,k-1))*(y1-ya)/(ya-y2)+uz(i+1,j,k-1)

                   bzxn(j,k)=r1*uu1/r2
                   bzyn(j,k)=   uv1
                   bzzn(j,k)=r1*uw1/r2
                enddo
             endif
          endif
       enddo
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Compute outflow
    call MPI_ALLREDUCE(inflow,outflow,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    outflow = outflow / nx / nz
    if (dims(1)==1) then
      j = xsize(2)
      do k = 1, xsize(3)
         do i = 1, xsize(1)
            byxn(i, k) = ux(i, j, k) - dt * outflow * (ux(i, j, k) - ux(i, j - 1, k))
            byyn(i, k) = uy(i, j, k) - dt * outflow * (uy(i, j, k) - uy(i, j - 1, k))
            byzn(i, k) = uz(i, j, k) - dt * outflow * (uz(i, j, k) - uz(i, j - 1, k))

            rho(i, j, k, 1) = rho(i, j, k, 1) &
                 - dt * outflow * (rho(i, j, k, 1) - rho(i, j - 1, k, 1))
         enddo
      enddo
    elseif (ny - (nym / dims(1)) == xstart(2)) then
       j = xsize(2) - 1
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, j, k) - dt * outflow * (ux(i, j, k) - ux(i, j - 1, k))
             byyn(i, k) = uy(i, j, k) - dt * outflow * (uy(i, j, k) - uy(i, j - 1, k))
             byzn(i, k) = uz(i, j, k) - dt * outflow * (uz(i, j, k) - uz(i, j - 1, k))

             rho(i, j, k, 1) = rho(i, j, k, 1) &
                  - dt * outflow * (rho(i, j, k, 1) - rho(i, j - 1, k, 1))
          enddo
       enddo
    endif

    return
  end subroutine boundary_conditions_jet

  !********************************************************************
  !
  subroutine jet_flrt (ux,constant)
    !
    !********************************************************************

    USE decomp_2d
    USE decomp_2d_poisson
    USE variables
    USE param
    USE var
    USE MPI

    implicit none

    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux
    real(mytype) :: constant

    integer :: j,i,k,code
    real(mytype) :: can,ut3,ut,ut4

    ut3=zero
    do k=1,ysize(3)
       do i=1,ysize(1)
          ut=zero
          do j=1,ny-1
             if (istret.eq.0) then
                ut=ut+dy*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             else
                ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-half*(ux(i,j+1,k)-ux(i,j,k)))
             endif
          enddo
          ut=ut/yly
          ut3=ut3+ut
       enddo
    enddo
    ut3=ut3/(real(nx*nz,mytype))

    call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can=-(constant-ut4)

    if (nrank==0) print *,nrank,'UT',ut4,can

    do k=1,ysize(3)
       do i=1,ysize(1)
          do j=2,ny-1
             ux(i,j,k)=ux(i,j,k)-can
          enddo
       enddo
    enddo

    return
  end subroutine jet_flrt
  !********************************************************************

  !############################################################################
  subroutine init_post(ep1)

    USE MPI

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1
    real(mytype) :: x, xprobes, yprobes, zprobes
    integer :: i,j,k,code
    character :: a

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post start'
#endif

    allocate(usum(ysize(2)),vsum(ysize(2)),wsum(ysize(2)))
    allocate(uusum(ysize(2)),uvsum(ysize(2)),uwsum(ysize(2)))
    allocate(vvsum(ysize(2)),vwsum(ysize(2)),wwsum(ysize(2)))
    usum=zero;vsum=zero;wsum=zero
    uusum=zero;uvsum=zero;uwsum=zero
    vvsum=zero;vwsum=zero;wwsum=zero
    ntimes1 = 0
    ntimes2 = 0
    nprobes  = 0

    !probes
    !WORK X-PENCILS
    open(10,file='probes.prm',status='unknown',form='formatted')
    read (10,*) nprobes
    read (10,*) a
    if (nprobes .gt. 0) then
       allocate(nxprobes(nprobes), nyprobes(nprobes), nzprobes(nprobes), rankprobes(nprobes))
       rankprobes(:)=0
       do i=1, nprobes
          read (10,*) xprobes, yprobes, zprobes
          !x
          if (nclx) then
             nxprobes(i)=int(xprobes/dx)
          else
             nxprobes(i)=int(xprobes/dx+1)
          end if
          !y
          if (ncly) then
             nyprobes(i)=int(yprobes/dy)
          else
             nyprobes(i)=int(yprobes/dy+1)
          end if
          !z
          if (nclz) then
             nzprobes(i)=int(zprobes/dz)
          else
             nzprobes(i)=int(zprobes/dz+1)
          end if
          if       (xstart(1) .le. nxprobes(i) .and. nxprobes(i) .le. xend(1)) then
             if    (xstart(2) .le. nyprobes(i) .and. nyprobes(i) .le. xend(2)) then
                if (xstart(3) .le. nzprobes(i) .and. nzprobes(i) .le. xend(3)) then
                   rankprobes(i)=1
                endif
             endif
          endif
       enddo
    endif
    close(10)

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init_post ok'
#endif

  end subroutine init_post
  !############################################################################
  subroutine postprocessing_jet(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    character(len=30) :: filename

    if (itime.ge.initstat) then
       !umean=ux1
       call fine_to_coarseS(1,ux1,tmean)
       umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)

       !vmean=uy1
       call fine_to_coarseS(1,uy1,tmean)
       vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

       !wmean=uz1
       call fine_to_coarseS(1,uz1,tmean)
       wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

       !uumean=ux1*ux1
       ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

       !vvmean=uy1*uy1
       ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

       !wwmean=uz1*uz1
       ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

       !uvmean=ux1*uy1
       ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

       !uwmean=ux1*uz1
       ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

       !vwmean=uy1*uz1
       ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
       call fine_to_coarseS(1,ta1,tmean)
       vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

       if (mod(itime,icheckpoint)==0) then
          if (nrank==0) print *,'===========================================================<<<<<'
          if (nrank==0) print *,'Writing stat file',itime
          write(filename,"('umean.dat',I7.7)") itime
          call decomp_2d_write_one(1,umean,filename,1)
          write(filename,"('vmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,vmean,filename,1)
          write(filename,"('wmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,wmean,filename,1)
          write(filename,"('uumean.dat',I7.7)") itime
          call decomp_2d_write_one(1,uumean,filename,1)
          write(filename,"('vvmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,vvmean,filename,1)
          write(filename,"('wwmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,wwmean,filename,1)
          write(filename,"('uvmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,uvmean,filename,1)
          write(filename,"('uwmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,uwmean,filename,1)
          write(filename,"('vwmean.dat',I7.7)") itime
          call decomp_2d_write_one(1,vwmean,filename,1)
          if (nrank==0) print *,'write stat done!'
          if (nrank==0) print *,'===========================================================<<<<<'
          if (nrank==0) then
             write(filename,"('umean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('vmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('wmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('uumean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('vvmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('wwmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('uvmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('uwmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
             write(filename,"('vwmean.dat',I7.7)") itime-icheckpoint
             call system ("rm " //filename)
          endif
       endif

    endif

    return
  end subroutine postprocessing_jet
  !############################################################################
  subroutine write_probes(ux1,uy1,uz1,phi1) !By Felipe Schuch

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux1, uy1, uz1
    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),numscalar) :: phi1

    integer :: i
    character(len=30) :: filename
    FS = 1+3+numscalar !Number of columns
    write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
    FS = FS*14+1  !Line width

    do i=1, nprobes
       if (rankprobes(i) .eq. 1) then
          write(filename,"('./out/probe',I4.4)") i
          open(67,file=trim(filename),status='unknown',form='formatted'&
               ,access='direct',recl=FS)
          write(67,fileformat,rec=itime) t,&                         !1
               ux1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !2
               uy1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !3
               uz1(nxprobes(i),nyprobes(i),nzprobes(i)),&            !4
               phi1(nxprobes(i),nyprobes(i),nzprobes(i),:),&         !numscalar
               NL                                                    !+1
          close(67)
       endif
    enddo

  end subroutine write_probes
  !############################################################################

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: momentum_forcing_jet
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies a fringe/sponge region at the outlet.
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE momentum_forcing_jet(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    INTEGER :: i, j, k
    REAL(mytype) :: y, yfringe
    REAL(mytype) :: f

    !! Set fringe height
    !! Fringe forcing will be applied for y > yfringe
    yfringe = 0.9_mytype * yly

    DO k = 1, xsize(3)
       DO j = 1, xsize(2)
          IF (istret.EQ.0) THEN
             y=REAL(j+xstart(2)-1-1,mytype)*dy
          ELSE
             y=yp(j+xstart(2)-1)
          ENDIF

          IF (y.GT.yfringe) THEN
             IF (y.LT.(yfringe + half * (yly - yfringe))) THEN
                f = SIN(((y - yfringe) / (yly - yfringe + 1.0e-16_mytype)) * (half * PI))
             ELSE
                f = one
             ENDIF

             DO i = 1, xsize(1)

                !! uy -> mean influx = outflow
                duy1(i, j, k, 1) = duy1(i, j, k, 1) + f * rho1(i, j, k, 1) * (outflow - uy1(i, j, k))

                !! ux,uz -> zero
                dux1(i, j, k, 1) = dux1(i, j, k, 1) + f * rho1(i, j, k, 1) * (zero - ux1(i, j, k))
                duz1(i, j, k, 1) = duz1(i, j, k, 1) + f * rho1(i, j, k, 1) * (zero - uz1(i, j, k))

             ENDDO
          ENDIF
       ENDDO
    ENDDO

  ENDSUBROUTINE momentum_forcing_jet

end module jet
