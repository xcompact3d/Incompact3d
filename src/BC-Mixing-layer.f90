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

module mixlayer

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: init_mixlayer!, boundary_conditions_mixlayer, postprocessing_mixlayer

contains

  subroutine init_mixlayer (rho1,ux1,uy1,uz1)

    use decomp_2d, only : mytype, xsize
    use param, only : u1, u2, dens1, dens2
    use param, only : half, one, two, four, eight, sixteen
    use param, only : ntime, nrhotime
    use dbg_schemes, only: sin_prec, cos_prec, exp_prec, tanh_prec, sqrt_prec
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
       u1 =  sqrt_prec(dens2 / dens1) / (sqrt_prec(dens2 / dens1) + one)
       u2 = -sqrt_prec(dens1 / dens2) / (one + sqrt_prec(dens1 / dens2))
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
                     + half * (u1 - u2) * tanh_prec(two * y)
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
                disturb_decay = 0.025_mytype * (u1 - u2) * exp_prec(-0.05_mytype * (y**2))
                u_disturb = disturb_decay * (sin_prec(eight * PI * x / xlx) &
                     + sin_prec(four * PI * x / xlx) / eight &
                     + sin_prec(two * PI * x / xlx) / sixteen)
                u_disturb = (0.05_mytype * y * xlx / PI) * u_disturb
                v_disturb = disturb_decay * (cos_prec(eight * PI * x / xlx) &
                     + cos_prec(four * PI * x / xlx) / eight &
                     + cos_prec(two * PI * x / xlx) / sixteen)

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
