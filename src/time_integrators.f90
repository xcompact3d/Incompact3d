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

module time_integrators

  implicit none

  private
  public :: int_time

contains

  subroutine intt(var1,dvar1,npaire,isc,forcing1)

    use MPI
    use param
    use variables
    use decomp_2d
    use ydiff_implicit, only : inttimp
#ifdef DEBG 
    use tools, only : avg3d
#endif

    implicit none

    !! INPUT / OUTPUT
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in), optional :: forcing1
    integer, intent(in), optional :: npaire, isc

    !! LOCAL
    integer :: is, code, ierror

#ifdef DEBG 
    real(mytype) avg_param
#endif

#ifdef DEBG
    avg_param = zero
    call avg3d (var1, avg_param)
    if (nrank == 0) write(*,*)'## SUB intt VAR var1 (start) AVG ', avg_param
    avg_param = zero
    call avg3d (dvar1(:,:,:,1), avg_param)
    if (nrank == 0) write(*,*)'## SUB intt VAR dvar1(1) (start) AVG ', avg_param
    avg_param = zero
    call avg3d (dvar1(:,:,:,2), avg_param)
    if (nrank == 0) write(*,*)'## SUB intt VAR dvar1(2) (start) AVG ', avg_param
#endif

    if (iimplicit.ge.1) then
       !>>> (semi)implicit Y diffusion

       if (present(isc)) then
          is = isc
       else
          is = 0
       endif
       if (present(npaire).and.present(forcing1)) then
          call inttimp(var1, dvar1, npaire=npaire, isc=is, forcing1=forcing1)
       else if (present(npaire)) then
          call inttimp(var1, dvar1, npaire=npaire, isc=is)
       else
          if (nrank  == 0) write(*,*) "Error in intt call."
          call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif

    elseif (itimescheme.eq.1) then
       !>>> Euler

       var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
    elseif(itimescheme.eq.2) then
       !>>> Adam-Bashforth second order (AB2)

       ! Do first time step with Euler
       if(itime.eq.1.and.irestart.eq.0) then
          var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
       else
          var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
       endif
       dvar1(:,:,:,2)=dvar1(:,:,:,1)
    elseif(itimescheme.eq.3) then
       !>>> Adams-Bashforth third order (AB3)

       ! Do first time step with Euler
       if(itime.eq.1.and.irestart.eq.0) then
          var1(:,:,:)=dt*dvar1(:,:,:,1)+var1(:,:,:)
       elseif(itime.eq.2.and.irestart.eq.0) then
          ! Do second time step with AB2
          var1(:,:,:)=onepfive*dt*dvar1(:,:,:,1)-half*dt*dvar1(:,:,:,2)+var1(:,:,:)
          dvar1(:,:,:,3)=dvar1(:,:,:,2)
       else
          ! Finally using AB3
          var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+cdt(itr)*dvar1(:,:,:,3)+var1(:,:,:)
          dvar1(:,:,:,3)=dvar1(:,:,:,2)
       endif
       dvar1(:,:,:,2)=dvar1(:,:,:,1)
    elseif(itimescheme.eq.4) then
       !>>> Adams-Bashforth fourth order (AB4)

       if (nrank==0) then
          write(*,*) "AB4 not implemented!"
          stop
       endif

       !if (itime.eq.1.and.ilit.eq.0) then
       !var(:,:,:)=gdt(itr)*hx(:,:,:)+var(:,:,:)
       !uy(:,:,:)=gdt(itr)*hy(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=gdt(itr)*hz(:,:,:)+uz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !elseif (itime.eq.2.and.ilit.eq.0) then
       !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+var(:,:,:)
       !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+uz(:,:,:)
       !gox(:,:,:)=gx(:,:,:)
       !goy(:,:,:)=gy(:,:,:)
       !goz(:,:,:)=gz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !elseif (itime.eq.3.and.ilit.eq.0) then
       !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+var(:,:,:)
       !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+uz(:,:,:)
       !gox(:,:,:)=gx(:,:,:)
       !goy(:,:,:)=gy(:,:,:)
       !goz(:,:,:)=gz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !else
       !var(:,:,:)=adt(itr)*hx(:,:,:)+bdt(itr)*gx(:,:,:)+cdt(itr)*gox(:,:,:)+ddt(itr)*gax(:,:,:)+var(:,:,:)
       !uy(:,:,:)=adt(itr)*hy(:,:,:)+bdt(itr)*gy(:,:,:)+cdt(itr)*goy(:,:,:)+ddt(itr)*gay(:,:,:)+uy(:,:,:)
       !uz(:,:,:)=adt(itr)*hz(:,:,:)+bdt(itr)*gz(:,:,:)+cdt(itr)*goz(:,:,:)+ddt(itr)*gaz(:,:,:)+uz(:,:,:)
       !gax(:,:,:)=gox(:,:,:)
       !gay(:,:,:)=goy(:,:,:)
       !gaz(:,:,:)=goz(:,:,:)
       !gox(:,:,:)=gx(:,:,:)
       !goy(:,:,:)=gy(:,:,:)
       !goz(:,:,:)=gz(:,:,:)
       !gx(:,:,:)=hx(:,:,:)
       !gy(:,:,:)=hy(:,:,:)
       !gz(:,:,:)=hz(:,:,:)
       !endif
       !>>> Runge-Kutta (low storage) RK3
    elseif(itimescheme.eq.5) then
       if(itr.eq.1) then
          var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
       else
          var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
       endif
       dvar1(:,:,:,2)=dvar1(:,:,:,1)
       !>>> Runge-Kutta (low storage) RK4
    elseif(itimescheme.eq.6) then

       if (nrank==0) then
          write(*,*) "RK4 not implemented!"
          STOP
       endif

    else

       if (nrank==0) then
          write(*,*) "Unrecognised itimescheme: ", itimescheme
          STOP
       endif

    endif

#ifdef DEBG
    avg_param = zero
    call avg3d (var1, avg_param)
    if (nrank == 0) write(*,*)'## SUB intt VAR var1 AVG ', avg_param
    avg_param = zero
    call avg3d (dvar1(:,:,:,1), avg_param)
    if (nrank == 0) write(*,*)'## SUB intt VAR dvar1 AVG ', avg_param
    if (nrank   ==  0) write(*,*)'# intt done'
#endif

    return

  end subroutine intt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time
  !! DESCRIPTION: 
  !!      INPUTS: 
  !!     OUTPUTS:
  !!       NOTES: 
  !!      AUTHOR:  
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE int_time(rho1, ux1, uy1, uz1, phi1, drho1, dux1, duy1, duz1, dphi1)

    use decomp_2d, only : mytype, xsize, nrank
    use param, only : zero, one
    use param, only : ntime, nrhotime, ilmn, iscalar, ilmn_solve_temp,itimescheme
    use param, only : iimplicit, sc_even
    use param, only : primary_species, massfrac
    use param, only : scalar_lbound, scalar_ubound
    use variables, only : numscalar,nu0nu
    use var, only : ta1, tb1
#ifdef DEBG 
    use tools, only : avg3d
#endif

    IMPLICIT NONE

    !! INPUT/OUTPUT
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1, dux1, duy1, duz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! OUTPUT
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

    !! LOCAL
    integer :: is, i, j, k
#ifdef DEBG
    real(mytype) avg_param
    if (nrank .eq. 0) write(*,*)'## Init int_time'
#endif

    call int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)
#ifdef DEBG
     avg_param = zero
     call avg3d (dux1, avg_param)
     if (nrank == 0) write(*,*)'## int_time dux1 ', avg_param
     avg_param = zero
     call avg3d (duy1, avg_param)
     if (nrank == 0) write(*,*)'## int_time duy1 ', avg_param
     avg_param = zero
     call avg3d (duz1, avg_param)
     if (nrank == 0) write(*,*)'## int_time duz1 ', avg_param
#endif

    IF (ilmn) THEN
       IF (ilmn_solve_temp) THEN
          CALL int_time_temperature(rho1, drho1, dphi1, phi1)
       ELSE
          CALL int_time_continuity(rho1, drho1)
       ENDIF
    ENDIF

    IF (iscalar.NE.0) THEN
       IF (ilmn.and.ilmn_solve_temp) THEN
          !! Compute temperature
          call calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))
       ENDIF

       DO is = 1, numscalar
          IF (is.NE.primary_species) THEN
             IF (iimplicit.ge.1) then
                if (sc_even(is)) then
                   k = 1
                else
                   k = 0
                endif
                CALL intt(phi1(:,:,:,is), dphi1(:,:,:,:,is), npaire=k, isc=is)
             ELSE
                CALL intt(phi1(:,:,:,is), dphi1(:,:,:,:,is))
             ENDIF

             DO k = 1, xsize(3)
                DO j = 1, xsize(2)
                   DO i = 1, xsize(1)
                      phi1(i,j,k,is) = max(phi1(i,j,k,is),scalar_lbound(is))
                      phi1(i,j,k,is) = min(phi1(i,j,k,is),scalar_ubound(is))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       IF (primary_species.GE.1) THEN
          phi1(:,:,:,primary_species) = one
          DO is = 1, numscalar
             IF ((is.NE.primary_species).AND.massfrac(is)) THEN
                phi1(:,:,:,primary_species) = phi1(:,:,:,primary_species) - phi1(:,:,:,is)
             ENDIF
          ENDDO

          DO k = 1, xsize(3)
             DO j = 1, xsize(2)
                DO i = 1, xsize(1)
                   phi1(i,j,k,primary_species) = max(phi1(i,j,k,primary_species),zero)
                   phi1(i,j,k,primary_species) = min(phi1(i,j,k,primary_species),one)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       IF (ilmn.and.ilmn_solve_temp) THEN
          !! Compute rho
          call calc_rho_eos(rho1(:,:,:,1), ta1, phi1, tb1, xsize(1), xsize(2), xsize(3))
       ENDIF
    ENDIF

#ifdef DEBG
    if (nrank .eq. 0) write(*,*)'## End  int_time'
#endif

  ENDSUBROUTINE int_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_momentum
  !! DESCRIPTION: Integrates the momentum equations in time by calling time
  !!              integrator.
  !!      INPUTS: dux1, duy1, duz1 - the RHS(s) of the momentum equations
  !!     OUTPUTS: ux1,   uy1,  uz1 - the intermediate momentum state.
  !!       NOTES: This is integrating the MOMENTUM in time (!= velocity)
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

    USE param
    USE variables
    USE var, ONLY: px1, py1, pz1
    USE decomp_2d

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1, duy1, duz1

    if (iimplicit.ge.1) then
       call intt(ux1, dux1, npaire=1, isc=0, forcing1=px1)
       call intt(uy1, duy1, npaire=0, isc=0, forcing1=py1)
       call intt(uz1, duz1, npaire=1, isc=0, forcing1=pz1)
    else
       call intt(ux1, dux1)
       call intt(uy1, duy1)
       call intt(uz1, duz1)
    endif

  endsubroutine int_time_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_continuity
  !! DESCRIPTION: Integrates the continuity (aka density transport) equation in
  !!              time
  !!      INPUTS: drho1 - the RHS(s) of the continuity equation.
  !!     OUTPUTS:  rho1 - the density at new time.
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_continuity(rho1, drho1)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    integer :: it, i, j, k
    real(mytype) :: rhomin, rhomax

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! First, update old density / store old transients depending on scheme
    if (itimescheme.lt.5) then
       !! Euler/AB - Store old density values
       do it = nrhotime, 2, -1
          rho1(:,:,:,it) = rho1(:,:,:,it-1)
       enddo
    elseif (itimescheme.eq.5) then
       !! RK3 - Stores old transients
       if (itr.eq.1) then
          do it = nrhotime, 2, -1
             rho1(:,:,:,it) = rho1(:,:,:,it-1)
          enddo
          rho1(:,:,:,2) = drho1(:,:,:,1)
       endif
    else
       if (nrank  == 0) then
          write(*,*) "int_time_continuity not implemented for itimescheme", itimescheme
          stop
       endif
    endif

    !! Now we can update current density
    call intt(rho1(:,:,:,1), drho1)

    !! Enforce boundedness on density
    if (ilmn_bound) then
       rhomin = min(dens1, dens2)
       rhomax = max(dens1, dens2)
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                rho1(i, j, k, 1) = max(rho1(i, j, k, 1), rhomin)
                rho1(i, j, k, 1) = min(rho1(i, j, k, 1), rhomax)
             enddo
          enddo
       enddo
    endif

  endsubroutine int_time_continuity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_temperature
  !! DESCRIPTION: Integrates the temperature equation in time
  !!      INPUTS: drho1 - the RHS(s) of the temperature equation.
  !!     OUTPUTS:  rho1 - the density at new time.
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_temperature(rho1, drho1, dphi1, phi1)

    USE param
    USE variables
    USE decomp_2d

    USE navier, ONLY : lmn_t_to_rho_trans
    USE var, ONLY : tc1, tb1

    implicit none

    integer :: it, i, j, k

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! First, update old density / store old transients depending on scheme
    if (itimescheme.lt.5) then
       !! Euler/AB - Store old density values
       do it = nrhotime, 2, -1
          rho1(:,:,:,it) = rho1(:,:,:,it-1)
       enddo
    elseif (itimescheme.eq.5) then
       !! RK3 - Stores old transients
       if (itr.eq.1) then
          do it = nrhotime, 2, -1
             rho1(:,:,:,it) = rho1(:,:,:,it-1)
          enddo

          !! Convert temperature transient to density transient and store it.
          call lmn_t_to_rho_trans(rho1(:,:,:,2), drho1(:,:,:,1), rho1(:,:,:,1), dphi1, phi1)
       endif
    else
       if (nrank==0) then
          write(*,*) "int_time_continuity not implemented for itimescheme", itimescheme
          stop
       endif
    endif

    !!-------------------------------------------------------------------
    !! XXX We are integrating the temperature equation - get temperature
    !!-------------------------------------------------------------------
    call calc_temp_eos(tc1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

    !! Now we can update current temperature
    call intt(tc1, drho1)

    !! Temperature >= 0
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             tc1(i,j,k) = max(tc1(i,j,k), zero)
          enddo
       enddo
    enddo

    !!-------------------------------------------------------------------
    !! XXX We are integrating the temperature equation - get back to rho
    !!-------------------------------------------------------------------
    call calc_rho_eos(rho1(:,:,:,1), tc1, phi1, tb1, xsize(1), xsize(2), xsize(3))

  endsubroutine int_time_temperature

end module time_integrators
