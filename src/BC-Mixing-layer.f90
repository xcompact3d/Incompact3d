!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!        FILE: BC-Mixing-layer
!!      AUTHOR: Paul Bartholomew
!! DESCRIPTION: A Low Mach, variable-density mixing layer with zero
!!              convective velocity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mixlayer

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: init_mixlayer!, boundary_conditions_mixlayer, postprocessing_mixlayer

contains

  subroutine init_mixlayer (rho1,ux1,uy1,uz1,phi1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : u1, u2, dens1, dens2
    USE param, ONLY : half, one, two, four, eight, sixteen
    USE param, ONLY : ntime, nrhotime
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i, j, k, is, it
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
       if (ilmn) then
          rhomin = MIN(dens1, dens2)
          rhomax = MAX(dens1, dens2)
          u1 = SQRT(dens2 / dens1) / (SQRT(dens2 / dens1) + one)
          u2 = -SQRT(dens1 / dens2) / (one + SQRT(dens1 / dens2))
       else
          rhomin = one
          rhomax = one
          u1 = half
          u2 = -half
       endif
       T1 = pressure0 / dens1
       T2 = pressure0 / dens2
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
                     + half * (u1 - u2) * TANH(two * y)
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
                disturb_decay = 0.025_mytype * (u1 - u2) * EXP(-0.05_mytype * (y**2))
                u_disturb = disturb_decay * (SIN(eight * PI * x / xlx) &
                     + SIN(four * PI * x / xlx) / eight &
                     + SIN(two * PI * x / xlx) / sixteen)
                u_disturb = (0.05_mytype * y * xlx / PI) * u_disturb
                v_disturb = disturb_decay * (COS(eight * PI * x / xlx) &
                     + COS(four * PI * x / xlx) / eight &
                     + COS(two * PI * x / xlx) / sixteen)

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
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_mixlayer

end module mixlayer
