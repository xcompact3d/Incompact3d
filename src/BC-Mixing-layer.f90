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
  
  subroutine init_mixlayer (rho1,ux1,uy1,uz1,drho1,dux1,duy1,duz1)

    USE decomp_2d, ONLY : mytype, xsize
    USE param, ONLY : u1, u2, dens1, dens2
    USE param, ONLY : half
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: rho1, drho1

    integer :: i, j, k, is
    real(mytype) :: x, y, z

    real(mytype) :: M, rspech, cp
    real(mytype) :: T1, T2, rhomin, rhomax
    real(mytype) :: disturb_decay, u_disturb, v_disturb

    integer, dimension (:), allocatable :: seed

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
       T1 = 1._mytype / dens1
       T2 = 1._mytype / dens2
       u1 = SQRT(dens2 / dens1) / (SQRT(dens2 / dens1) + 1._mytype)
       u2 = -SQRT(dens1 / dens2) / (1._mytype + SQRT(dens1 / dens2))
       M = 0.2_mytype
       rspech = 1.4_mytype
       cp = (1._mytype / (T2 * (rspech - 1._mytype))) * ((u1 - u2) / M)**2

       do k=1,xsize(3)
          do j=1,xsize(2)
             y=real((j+xstart(2)-1-1),mytype)*dy - half * yly
             do i=1,xsize(1)
                x=real(i-1,mytype)*dx

                !! Set mean field
                ux1(i, j, k) = ux1(i, j, k) + (u1 + u2) / 2._mytype &
                     + (u1 - u2) * TANH(2._mytype * y) / 2._mytype
                uy1(i, j, k) = zero
                uz1(i, j, k) = zero
                
                rho1(i, j, k, 1) = (1._mytype / (2._mytype * cp)) &
                     * (ux1(i, j, k) * (u1 + u2) - ux1(i, j, k)**2 - u1 * u2) &
                     + ux1(i, j, k) * (T1 - T2) / (u1 - u2) &
                     + (T2 * u1 - T1 * u2) / (u1 - u2)
                rho1(i, j, k, 1) = 1._mytype / rho1(i, j, k, 1)
                rho1(i, j, k, 1) = MAX(rho1(i, j, k, 1), rhomin)
                rho1(i, j, k, 1) = MIN(rho1(i, j, k, 1), rhomax)

                ! Calculate disturbance field (as given in Fortune2004)
                ! NB x and y are swapped relative to Fortune2004
                disturb_decay = 0.025 * (u1 - u2) * EXP(-0.05_mytype * (y**2))
                u_disturb = disturb_decay * (SIN(8._mytype * PI * x / xlx) &
                     + SIN(4._mytype * PI * x / xlx) / 8._mytype &
                     + SIN(2._mytype * PI * x / xlx) / 16._mytype)
                u_disturb = (0.05_mytype * y * xlx / PI) * u_disturb
                v_disturb = disturb_decay * (COS(8._mytype * PI * x / xlx) &
                     + COS(4._mytype * PI * x / xlx) / 8._mytype &
                     + COS(2._mytype * PI * x / xlx) / 16._mytype)

                ux1(i, j, k) = ux1(i, j, k) + u_disturb
                uy1(i, j, k) = uy1(i, j, k) + v_disturb
             enddo
          enddo
       enddo

    endif

    dux1(:,:,:,1)=ux1(:,:,:)
    duy1(:,:,:,1)=uy1(:,:,:)
    duz1(:,:,:,1)=uz1(:,:,:)
    do is = 2, ntime
       dux1(:,:,:,is)=dux1(:,:,:,is - 1)
       duy1(:,:,:,is)=duy1(:,:,:,is - 1)
       duz1(:,:,:,is)=duz1(:,:,:,is - 1)
    enddo
    
    drho1(:,:,:,1)=rho1(:,:,:,1)
    do is = 2, ntime
       rho1(:,:,:,is)=rho1(:,:,:,is - 1)
       drho1(:,:,:,is)=drho1(:,:,:,is - 1)
    enddo

#ifdef DEBG
    if (nrank .eq. 0) print *,'# init end ok'
#endif

    return
  end subroutine init_mixlayer
  
end module mixlayer
