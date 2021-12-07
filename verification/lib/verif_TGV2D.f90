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

module verif_tgv2D

  USE decomp_2d
  USE variables
  USE param
  USE tools, only : error_l1_l2_linf

  IMPLICIT NONE

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: error_tgv2D

contains

  !############################################################################
  subroutine error_tgv2D(ux1, uy1, phi1, errx,   erry,   errs, &
                                         errxd1, erryd1, errsd1, &
                                         errxd2, erryd2, errsd2)

    use decomp_2d
    use MPI
    use param, only : one, two, xnu, ifirst, itime
    use variables, only : numscalar, sc
    use dbg_schemes, only: sin_prec, cos_prec

    implicit none

    ! Arguments
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    ! Analytical error
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: errx, erry
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errs
    ! Time-discrete error
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: errxd1, erryd1
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errsd1
    ! Space-time discrete error
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: errxd2, erryd2
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: errsd2

    ! Local variables
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tmp
    integer :: i,j,k,l,it
    real(mytype) :: x, y, xl1, xl2, xlinf, yl1, yl2, ylinf
    real(mytype) :: solx0, soly0, sols0
    real(mytype) :: solxt, solyt, solst, solxd, solyd, solsd, solxdd, solydd, solsdd
    real(mytype) :: xdamping(3), ydamping(3), sdamping(3,numscalar)

    ! Compute the errors
    call compute_tgv2D_errors(xdamping, ydamping, sdamping)

    ! Compute solutions and errors
    do k = 1,xsize(3)
       do j = 1,xsize(2)
          y = real(j+xstart(2)-1-1,mytype)*dy
          do i = 1,xsize(1)
            x = real(i+xstart(1)-1-1,mytype)*dx
            ! Initial solution
            solx0 = sin_prec(x) * cos_prec(y)
            soly0 = - cos_prec(x) * sin_prec(y)
            ! Analytical solution
            solxt = solx0 * xdamping(1)
            solyt = soly0 * ydamping(1)
            ! Time discrete solution
            solxd = solx0 * xdamping(2)
            solyd = soly0 * ydamping(2)
            ! Space-time discrete solution
            solxdd = solx0 * xdamping(3)
            solydd = soly0 * ydamping(3)
            ! Errors
            errx(i,j,k) = ux1(i,j,k) - solxt
            erry(i,j,k) = uy1(i,j,k) - solyt
            errxd1(i,j,k) = ux1(i,j,k) - solxd
            erryd1(i,j,k) = uy1(i,j,k) - solyd
            errxd2(i,j,k) = ux1(i,j,k) - solxdd
            erryd2(i,j,k) = uy1(i,j,k) - solydd
          enddo
       enddo
    enddo
    if (iscalar==1) then
       do l = 1, numscalar
          do k = 1,xsize(3)
             do j = 1,xsize(2)
                y = real(j+xstart(2)-1-1,mytype)*dy
                do i = 1,xsize(1)
                  x = real(i+xstart(1)-1-1,mytype)*dx
                  ! Initial solution
                  sols0 = sin_prec(x) * sin_prec(y)
                  ! Analytical solution
                  solst = sols0 * sdamping(1,l)
                  ! Time discrete solution
                  solsd = sols0 * sdamping(2,l)
                  ! Space-time discrete solution
                  solsdd = sols0 * sdamping(3,l)
                  ! Errors
                  errs(i,j,k,l) = phi1(i,j,k,l) - solst
                  errsd1(i,j,k,l) = phi1(i,j,k,l) - solsd
                  errsd2(i,j,k,l) = phi1(i,j,k,l) - solsdd
                enddo
             enddo
          enddo
       enddo
    endif

  end subroutine error_tgv2D

  ! Compute the damping factors
  subroutine compute_tgv2D_errors(xdamping, ydamping, sdamping)

    use decomp_2d
    use param, only : one, two, xnu, ifirst, itime, itimescheme, iimplicit
    use variables, only : numscalar, sc
    use dbg_schemes, only: exp_prec

    implicit none

    real(mytype), intent(out) :: xdamping(3), ydamping(3), sdamping(3, numscalar)

    integer :: it, l
    real(mytype) :: ktgv, k2tgv, coef(numscalar+1)

    ! Compute modified wavenumber
    ktgv = one
    call compute_k2(ktgv, k2tgv)

    ! Init
    xdamping(:) = one
    ydamping(:) = one
    sdamping(:,:) = one

    ! Compute analytical damping
    coef(1) = exp_prec(-two*dt*xnu)
    do it = ifirst, itime
       xdamping(1) = xdamping(1) * coef(1)
       ydamping(1) = ydamping(1) * coef(1)
       do l = 1, numscalar
         sdamping(1,l) = sdamping(1,l) * exp_prec(-two*dt*xnu/sc(l))
       enddo
    enddo

    ! Compute time and space-time discrete damping
    !
    ! Explicit Euler
    if (itimescheme == 1) then

      ! Time discrete errors
      if (iimplicit == 0) then
        coef(1) = (one - two*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = one - two*dt*xnu/sc(l)
        enddo
      else if (iimplicit == 1) then
        coef(1) = (one - dt*xnu) / (one + dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - dt*xnu/sc(l)) / (one + dt*xnu/sc(l))
        enddo
      else if (iimplicit == 2) then
        coef(1) = (one - onepfive*dt*xnu) / (one + half*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - onepfive*dt*xnu/sc(l)) / (one + half*dt*xnu/sc(l))
        enddo
      endif
      do it = ifirst, itime
        xdamping(2) = xdamping(2) * coef(1)
        ydamping(2) = ydamping(2) * coef(1)
        do l = 1, numscalar
          sdamping(2,l) = sdamping(2,l) * coef(1+l)
        enddo
      enddo

      ! Space-time discrete errors
      if (iimplicit == 0) then
        coef(1) = (one - two*k2tgv*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = one - two*k2tgv*dt*xnu/sc(l)
        enddo
      else if (iimplicit == 1) then
        coef(1) = (one - k2tgv*dt*xnu) / (one + k2tgv*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - k2tgv*dt*xnu/sc(l)) / (one + k2tgv*dt*xnu/sc(l))
        enddo
      else if (iimplicit == 2) then
        coef(1) = (one - onepfive*k2tgv*dt*xnu) / (one + half*k2tgv*dt*xnu)
        do l = 1, numscalar
          coef(1+l) = (one - onepfive*k2tgv*dt*xnu/sc(l)) / (one + half*k2tgv*dt*xnu/sc(l))
        enddo
      endif
      do it = ifirst, itime
        xdamping(3) = xdamping(3) * coef(1)
        ydamping(3) = ydamping(3) * coef(1)
        do l = 1, numscalar
          sdamping(3,l) = sdamping(3,l) * coef(1+l)
        enddo
      enddo

    else

      if (nrank==0) write(*,*)  "TGV2D: No discrete error implemented for this time scheme."
      xdamping(2:) = xdamping(1)
      ydamping(2:) = ydamping(1)
      do l = 1, numscalar
        sdamping(2:,l) = sdamping(1,l)
      enddo

    endif

  end subroutine compute_tgv2D_errors

  ! Compute the modified wavenumber for the second derivative
  ! Warning : we use the X momentum wavenumber for Y momentum and for the scalars
  subroutine compute_k2(kin, k2out)

    use decomp_2d, only : mytype
    use param
    use derivx, only : alsaix, asix, bsix, csix, dsix
    use dbg_schemes, only: cos_prec

    implicit none

    real(mytype), intent(in) :: kin
    real(mytype), intent(out) :: k2out

    if (kin < zero .or. kin > pi/min(dx,dy)) then
      if (nrank==0) then
        write(*,*) "TGV2D: Warning, incorrect wavenumber provided."
      endif
    endif

    k2out = asix * two * (one - cos_prec(kin*dx)) &
          + four * bsix * half * (one - cos_prec(two*kin*dx)) &
          + nine * csix * (two / nine) * (one - cos_prec(three*kin*dx)) &
          + sixteen * dsix * (one / eight) * (one - cos_prec(four*kin*dx))
    k2out = k2out / (one + two * alsaix * cos_prec(kin*dx))

  end subroutine compute_k2

end module verif_tgv2D
