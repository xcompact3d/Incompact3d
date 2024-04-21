!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module mixlayer

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: init_mixlayer!, boundary_conditions_mixlayer, postprocessing_mixlayer

contains

  subroutine init_mixlayer (rho1,ux1,uy1,uz1)

    use param, only : u1, u2, dens1, dens2
    use param, only : half, one, two, four, eight, sixteen
    use param, only : ntime, nrhotime
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    integer :: i, j, k, is
    real(mytype) :: x, y, z

    real(mytype) :: M, rspech, heatcap
    real(mytype) :: T1, T2, rhomin, rhomax
    real(mytype) :: disturb_decay, u_disturb, v_disturb

    if (iin.eq.0) then !empty domain
       if (nrank==0) then
          write(*,*) "Empty initial domain!"
       endif

       ux1=zero; uy1=zero; uz1=zero
    endif

    if (iin.eq.1) then !generation of a random noise
       if (nrank==0) then
          write(*,*) "Filled initial domain!"
       endif

       ux1=zero; uy1=zero; uz1=zero

       !! Compute flow for zero convective velocity
       rhomin = min(dens1, dens2)
       rhomax = max(dens1, dens2)
       T1 = pressure0 / dens1
       T2 = pressure0 / dens2
       u1 =  sqrt(dens2 / dens1) / (sqrt(dens2 / dens1) + one)
       u2 = -sqrt(dens1 / dens2) / (one + sqrt(dens1 / dens2))
       M = 0.2_mytype
       rspech = 1.4_mytype
       heatcap = (one / (T2 * (rspech - one))) * ((u1 - u2) / M)**2

       do k=1,xsize(3)
          do j=1,xsize(2)
             y=real((j+xstart(2)-2),mytype)*dy - half * yly
             do i=1,xsize(1)
                x=real(i+xstart(1)-2,mytype)*dx

                !! Set mean field
                ux1(i, j, k) = ux1(i, j, k) + half * (u1 + u2) &
                     + half * (u1 - u2) * tanh(two * y)
                uy1(i, j, k) = zero
                uz1(i, j, k) = zero

                rho1(i, j, k, 1) = (one / (two * heatcap)) &
                     * (ux1(i, j, k) * (u1 + u2) - ux1(i, j, k)**2 - u1 * u2) &
                     + ux1(i, j, k) * (T1 - T2) / (u1 - u2) &
                     + (T2 * u1 - T1 * u2) / (u1 - u2)
                rho1(i, j, k, 1) = one / rho1(i, j, k, 1)
                rho1(i, j, k, 1) = MAX(rho1(i, j, k, 1), rhomin)
                rho1(i, j, k, 1) = MIN(rho1(i, j, k, 1), rhomax)

                ! Calculate disturbance field (as given in Fortune2004)
                disturb_decay = 0.025_mytype * (u1 - u2) * exp(-0.05_mytype * (y**2))
                u_disturb = disturb_decay * (sin(eight * PI * x / xlx) &
                     + sin(four * PI * x / xlx) / eight &
                     + sin(two * PI * x / xlx) / sixteen)
                u_disturb = (0.05_mytype * y * xlx / PI) * u_disturb
                v_disturb = disturb_decay * (cos(eight * PI * x / xlx) &
                     + cos(four * PI * x / xlx) / eight &
                     + cos(two * PI * x / xlx) / sixteen)

                ux1(i, j, k) = ux1(i, j, k) + u_disturb
                uy1(i, j, k) = uy1(i, j, k) + v_disturb
             enddo
          enddo
       enddo

       if (.not.ilmn) then
          rho1(:,:,:,:) = one
       endif

    endif

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_mixlayer

end module mixlayer
