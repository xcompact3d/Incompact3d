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

MODULE case

  USE param
  USE decomp_2d
  USE variables

  USE user_sim
  USE tgv
  USE cyl
  USE dbg_schemes
  USE channel
  USE mixlayer
  USE lockexch
  USE tbl

  USE var, ONLY : nzmsize

  IMPLICIT NONE

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init, boundary_conditions, postprocess_case, momentum_forcing, set_fluid_properties

CONTAINS

  SUBROUTINE init (rho1, ux1, uy1, uz1, ep1, phi1, drho1, dux1, duy1, duz1, dphi1, &
       pp3, px1, py1, pz1)

    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1,drho1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    REAL(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1

    INTEGER :: it, is

    !! Zero out the pressure field
    pp3(:,:,:,1) = zero
    px1(:,:,:) = zero
    py1(:,:,:) = zero
    pz1(:,:,:) = zero

    !! Default density and pressure0 to one
    pressure0 = one
    rho1(:,:,:,:) = one

    IF (itype.EQ.itype_user) THEN

       CALL init_user (ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_lockexch) THEN

       CALL init_lockexch(rho1, ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_tgv) THEN

       CALL init_tgv (ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_channel) THEN

       CALL init_channel (ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_cyl) THEN

       CALL init_cyl (ux1, uy1, uz1, phi1)

    ELSEIF (itype.EQ.itype_dbg) THEN

       CALL init_dbg (ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_tbl) THEN

       CALL init_tbl (ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_mixlayer) THEN

       CALL init_mixlayer (rho1, ux1, uy1, uz1, phi1)

    ELSE

       if (nrank.eq.0) then
          print *, "ERROR: Unknown itype: ", itype
          STOP
       endif

    ENDIF

    !! Setup old arrays
    do it = 1, ntime
       drho1(:,:,:,it) = rho1(:,:,:,1)
       dux1(:,:,:,it)=ux1(:,:,:)
       duy1(:,:,:,it)=uy1(:,:,:)
       duz1(:,:,:,it)=uz1(:,:,:)
    enddo

    do it = 2, nrhotime
       rho1(:,:,:,it) = rho1(:,:,:,1)
    enddo

    do is = 1, numscalar
       do it = 1, ntime
          dphi1(:,:,:,it,is) = phi1(:,:,:,is)
       enddo
    enddo

  END SUBROUTINE init

  SUBROUTINE boundary_conditions (rho,ux,uy,uz,phi,ep)

    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),nrhotime) :: rho

    IF (itype.EQ.itype_user) THEN

       CALL boundary_conditions_user (ux,uy,uz,phi,ep)

    ELSEIF (itype.EQ.itype_lockexch) THEN

       CALL boundary_conditions_lockexch(rho, phi)

    ELSEIF (itype.EQ.itype_tgv) THEN

       CALL boundary_conditions_tgv (ux, uy, uz, phi)

    ELSEIF (itype.EQ.itype_channel) THEN

       CALL boundary_conditions_channel (ux, uy, uz, phi)

    ELSEIF (itype.EQ.itype_cyl) THEN

       CALL boundary_conditions_cyl (ux, uy, uz, phi)

    ELSEIF (itype.EQ.itype_dbg) THEN

       CALL boundary_conditions_dbg (ux, uy, uz, phi)

    ELSEIF (itype.EQ.itype_tbl) THEN

       CALL boundary_conditions_tbl (ux, uy, uz, phi)

    ENDIF

  END SUBROUTINE boundary_conditions

  SUBROUTINE postprocess_case(rho,ux,uy,uz,pp,phi,ep)

    USE forces
    USE var, ONLY : nzmsize
    USE param, ONLY : npress

    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ep
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp

    IF (itype.EQ.itype_user) THEN

       CALL postprocess_user (ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_lockexch) THEN

       CALL postprocess_lockexch(rho, ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_tgv) THEN

       CALL postprocess_tgv (ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_channel) THEN

       CALL postprocess_channel (ux, uy, uz, pp, phi, ep)

    ELSEIF (itype.EQ.itype_cyl) THEN

       CALL postprocess_cyl (ux, uy, uz, ep)

    ELSEIF (itype.EQ.itype_dbg) THEN

       CALL postprocess_dbg (ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_tbl) THEN

       CALL postprocess_tbl (ux, uy, uz, ep)

    ENDIF

    if (iforces.eq.1) then
       call force(ux,uy,ep)
       call restart_forces(1)
    endif

  END SUBROUTINE postprocess_case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Calls case-specific forcing functions for the
  !!              momentum equations.
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE momentum_forcing(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

    IMPLICIT NONE

    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    REAL(mytype), INTENT(IN), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    IF (itype.EQ.itype_channel) THEN

       CALL momentum_forcing_channel(dux1, duy1, ux1, uy1)

    ENDIF

  ENDSUBROUTINE momentum_forcing

  subroutine set_fluid_properties(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    if (itype.eq.itype_lockexch) then

       call set_fluid_properties_lockexch(rho1, mu1)

    endif

  endsubroutine set_fluid_properties

END MODULE case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
