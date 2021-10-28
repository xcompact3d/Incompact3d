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

  REAL(mytype) :: outflow

  LOGICAL :: initialising

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_jet, boundary_conditions_jet, postprocess_jet, momentum_forcing_jet

contains

  subroutine init_jet (rho1,ux1,uy1,uz1,ep1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE var, ONLY : rho2, phi2, ta2, tb2, tc2, td2
    use dbg_schemes, only: exp_prec

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y
    integer :: i, j, k, ii, is, it, code
    real(mytype) :: um

    integer, dimension (:), allocatable :: seed

    ux1=zero
    uy1=zero
    uz1=zero
    if (iin /= 0) then
       call system_clock(count=code)
       if (iin == 2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
    endif

    !modulation of the random noise + initial velocity profile
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)
          um=exp_prec(-y**2)
          um = one
          do i=1,xsize(1)
             ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)
             uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
             uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
          enddo
       enddo
    enddo

    initialising = .true.
    call boundary_conditions_jet (rho1,ux1,uy1,uz1,phi1)
    initialising = .false.

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
    call transpose_x_to_y(rho1(:,:,:,1), rho2)
    do is = 1, numscalar
       call transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))
    enddo
    do k=1,ysize(3)
       do j=2,ysize(2)
          do i=1,ysize(1)
             ta2(i,j,k)=ta2(i,j,k)+ta2(i,1,k)
             tb2(i,j,k)=tb2(i,j,k)+tb2(i,1,k)
             tc2(i,j,k)=tc2(i,j,k)+tc2(i,1,k)

             rho2(i,j,k)=rho2(i,1,k)
             phi2(i,j,k,:) = phi2(i,1,k,:)
          enddo
       enddo
    enddo
    call transpose_y_to_x(ta2, ux1)
    call transpose_y_to_x(tb2, uy1)
    call transpose_y_to_x(tc2, uz1)
    call transpose_y_to_x(rho2, rho1(:,:,:,1))
    do is = 1, numscalar
       call transpose_y_to_x(phi2(:,:,:,is), phi1(:,:,:,is))
    enddo

#ifdef DEBG
    if (nrank == 0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_jet
  !********************************************************************
  subroutine boundary_conditions_jet (rho,ux,uy,uz,phi)

    USE MPI
    USE param
    USE variables
    USE decomp_2d
    USE var, only : ta1, tb1
    use dbg_schemes, only: sin_prec, tanh_prec, sqrt_prec

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho

    integer :: i, j, k, is
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

    if (t < one) then
       timeswitch = sin_prec(t * pi * half)
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
             r = sqrt_prec(x**2 + z**2) + 1.e-16_mytype !! Can't divide by zero!

             !! Set mean inflow
             byx1(i, k) = zero
             byy1(i, k) = (u1 - u2) * half &
                  * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D))) + u2
             byz1(i, k) = zero

             if (ilmn) then
                if (.not.ilmn_solve_temp) then
                   rho(i, 1, k, 1) = (dens1 - dens2) * half &
                        * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D))) + dens2
                else
                   ta1(i,1,k) = one !! Could set variable temperature inlet here
                endif
             else
                rho(i, 1, k, 1) = one
             endif

             if (iscalar/=0) then
                do is = 1, numscalar
                   if (.not.massfrac(is)) then
                      phi(i, 1, k, is) = half * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D)))
                   else if (is/=primary_species) then
                      phi(i, 1, k, is) = one - half * (one + tanh_prec((12.5_mytype / four) * ((D / two) / r - two * r / D)))
                   endif
                enddo

                if (primary_species.gt.0) then
                   phi(i, 1, k, primary_species) = one

                   do is = 1, numscalar
                      if (massfrac(is).and.(is.ne.primary_species)) then
                         phi(i, 1, k, primary_species) = phi(i, 1, k, primary_species) - phi(i, 1, k, is)
                      endif
                   enddo
                endif
             endif

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

       if (ilmn.and.ilmn_solve_temp) then
          !! Need to compute rho (on boundary)
          CALL calc_rho_eos(rho(:,1,:,1), ta1(:,1,:), phi(:,1,:,:), tb1(:,1,:), xsize(1), 1, xsize(3))
       endif
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

    if (initialising) then !! we can stop here
       return
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Apply lateral boundary conditions
    !! XXX Assume purely radial flow
    !! XXX In X-pencils

    xc=half*xlx
    zc=half*xlx

    !! X-BC
    IF (nclx1 == 2) THEN
       if(xstart(1) == 1)then!
          x = -xc
          i = 1
          do k=1,xsize(3)
             z=real(k + xstart(3) - 2,mytype)*dz-zc
             x2=x
             y2=z
             x1=x+dx
             y1=y2*x1/x2
             r1=sqrt_prec(x1**2+y1**2)
             r2=sqrt_prec(x2**2+y2**2)
             if(r1.gt.r2) write(*,*)'##BUG CL in JET'
             if(k==1)then!cas premier point
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
                if (z.gt.0.)then
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

          if (iscalar.ne.0) then
             do is = 1, numscalar
                if (is.ne.primary_species) then
                   phi(i,:,:,is) = one
                else
                   phi(i,:,:,is) = zero
                endif
             enddo
          endif

          if (ilmn) then
             if (.not.ilmn_solve_temp) then
                rho(i,:,:,1) = dens1
             else
                ta1(i,:,:) = one

                !! Need to compute rho (on boundary)
                CALL calc_rho_eos(rho(i,:,:,1), ta1(i,:,:), phi(i,:,:,:), tb1(i,:,:), 1, xsize(2), xsize(3))
             endif
          endif
       endif
    endif

    if (nclxn.eq.2) then
       if(xend(1).eq.nx)then
          x=xc
          i = xsize(1)
          do k=1,xsize(3)
             z=real(k + xstart(3) - 2,mytype)*dz-zc
             x2=x
             y2=z
             x1=x-dx
             y1=y2*x1/x2
             r1=sqrt_prec(x1**2+y1**2)
             r2=sqrt_prec(x2**2+y2**2)
             if(r1 > r2) write(*,*) '##Bug CL Jet 2'
             if (k == 1) then!cas premier point
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

          if (iscalar.ne.0) then
             do is = 1, numscalar
                if (is.ne.primary_species) then
                   phi(i,:,:,is) = one
                else
                   phi(i,:,:,is) = zero
                endif
             enddo
          endif

          if (ilmn) then
             if (.not.ilmn_solve_temp) then
                rho(i,:,:,1) = dens1
             else
                ta1(i,:,:) = one

                !! Need to compute rho (on boundary)
                call calc_rho_eos(rho(i,:,:,1), ta1(i,:,:), phi(i,:,:,:), tb1(i,:,:), 1, xsize(2), xsize(3))
             endif
          endif
       endif
    endif

    !! Z-BC
    if ((nclz1 == 2).and.(xstart(3) == 1)) then
       k = 1
       z = -zc
       do i=1,xsize(1)
          x=real(i + xstart(1) - 2, mytype) * dx - xc
          x2=z
          y2=x
          x1=z+dz
          y1=y2*x1/x2
          r1=sqrt_prec(x1**2+y1**2)
          r2=sqrt_prec(x2**2+y2**2)
          if(r1 > r2) write(*,*) 'bug CL'
          if(i == 1)then!cas premier point
             do j=1,xsize(2)

                bzx1(i,j)=r1*ux(i+1,j,k + 1)/r2
                bzy1(i,j)=   uy(i+1,j,k + 1)
                bzz1(i,j)=r1*uz(i+1,j,k + 1)/r2
             enddo
          elseif(i == (nx-1)/2+1)then!cas point du milieu
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

       if (iscalar.ne.0) then
          do is = 1, numscalar
             if (is.ne.primary_species) then
                phi(:,:,k,is) = one
             else
                phi(:,:,k,is) = zero
             endif
          enddo
       endif

       if (ilmn) then
          if (.not.ilmn_solve_temp) then
             rho(:,:,k,1) = dens1
          else
             ta1(:,:,k) = one

             !! Need to compute rho (on boundary)
             CALL calc_rho_eos(rho(:,:,k,1), ta1(:,:,k), phi(:,:,k,:), tb1(:,:,k), xsize(1), xsize(2), 1)
          endif
       endif
    endif

    if ((nclzn == 2).and.(xend(3) == nz)) then
       z=zc
       k = xsize(3)
       do i=1,xsize(1)
          x=real(i + xstart(1) - 2, mytype)*dx-xc
          x2=z
          y2=x
          x1=z-dz
          y1=y2*x1/x2
          r1=sqrt_prec(x1**2+y1**2)
          r2=sqrt_prec(x2**2+y2**2)
          if(r1 > r2) write(*,*) 'bug CL'
          if(i == 1)then!cas premier point
             do j=1,xsize(2)

                bzxn(i,j)=r1*ux(i+1,j,k-1)/r2
                bzyn(i,j)=   uy(i+1,j,k-1)
                bzzn(i,j)=r1*uz(i+1,j,k-1)/r2
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

       if (iscalar.ne.0) then
          do is = 1, numscalar
             if (is.ne.primary_species) then
                phi(:,:,k,is) = one
             else
                phi(:,:,k,is) = zero
             endif
          enddo
       endif

       if (ilmn) then
          if (.not.ilmn_solve_temp) then
             rho(:,:,k,1) = dens1
          else
             ta1(:,:,k) = one

             !! Need to compute rho (on boundary)
             CALL calc_rho_eos(rho(:,:,k,1), ta1(:,:,k), phi(:,:,k,:), tb1(:,:,k), xsize(1), xsize(2), 1)
          endif
       endif
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Compute outflow
    call MPI_ALLREDUCE(inflow,outflow,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    outflow = outflow / nx / nz
    if (xend(2).eq.ny) then
       j = xsize(2)
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, j, k) - dt * outflow * (ux(i, j, k) - ux(i, j - 1, k))
             byyn(i, k) = uy(i, j, k) - dt * outflow * (uy(i, j, k) - uy(i, j - 1, k))
             byzn(i, k) = uz(i, j, k) - dt * outflow * (uz(i, j, k) - uz(i, j - 1, k))
          enddo
       enddo

       if (iscalar.ne.0) then
          phi(:,j,:,:) = phi(:,j,:,:) - dt * outflow * (phi(:,j,:,:) - phi(:,j - 1,:,:))

          if (primary_species.gt.0) then
             phi(:,j,:,primary_species) = one
             do is = 1, numscalar
                if (massfrac(is).and.(is.ne.primary_species)) then
                   phi(:,j,:,primary_species) = phi(:,j,:,primary_species) - phi(:,j,:,is)
                endif
             enddo
          endif
       endif

       if (ilmn) then
          if (.not.ilmn_solve_temp) then
             rho(:, j, :, 1) = rho(:, j, :, 1) &
                  - dt * outflow * (rho(:, j, :, 1) - rho(:, j - 1, :, 1))
          else
             !! Compute temperature at j-1:j to form advection equation
             CALL calc_temp_eos(ta1(:,j-1:j,:), rho(:,j-1:j,:,1), phi(:,j-1:j,:,:), tb1(:,j-1:j,:), xsize(1), 2, xsize(3))

             ta1(:,j,:) = ta1(:,j,:) - dt * outflow * (ta1(:,j,:) - ta1(:,j-1,:))

             !! Need to compute rho (on boundary)
             CALL calc_rho_eos(rho(:,j,:,1), ta1(:,j,:), phi(:,j,:,:), tb1(:,j,:), xsize(1), 1, xsize(3))
          endif
       endif
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

  subroutine init_post(ep1)

    real(mytype),intent(in),dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ep1

  end subroutine init_post

  !############################################################################
  subroutine postprocess_jet(ux1,uy1,uz1,phi1,ep1) !By Felipe Schuch

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

    return
  end subroutine postprocess_jet
  !############################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: momentum_forcing_jet
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies a fringe/sponge region at the outlet.
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine momentum_forcing_jet(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

    use dbg_schemes, only: sin_prec

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    integer :: i, j, k
    real(mytype) :: y, yfringe
    real(mytype) :: f

    !! Set fringe height
    !! Fringe forcing will be applied for y > yfringe
    yfringe = 0.9_mytype * yly

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          if (istret.eq.0) then
             y=real(j+xstart(2)-1-1,mytype)*dy
          else
             y=yp(j+xstart(2)-1)
          endif

          if (y.gt.yfringe) then
             if (y.lt.(yfringe + half * (yly - yfringe))) then
                f = sin_prec(((y - yfringe) / (yly - yfringe + 1.0e-16_mytype)) * (half * pi))
             else
                f = one
             endif

             do i = 1, xsize(1)

                !! uy -> mean influx = outflow
                duy1(i, j, k, 1) = duy1(i, j, k, 1) + f * rho1(i, j, k, 1) * (outflow - uy1(i, j, k))

                !! ux,uz -> zero
                dux1(i, j, k, 1) = dux1(i, j, k, 1) + f * rho1(i, j, k, 1) * (zero - ux1(i, j, k))
                duz1(i, j, k, 1) = duz1(i, j, k, 1) + f * rho1(i, j, k, 1) * (zero - uz1(i, j, k))

             enddo
          endif
       enddo
    enddo

  endsubroutine momentum_forcing_jet

end module jet
