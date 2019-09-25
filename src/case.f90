!! -*- mode: F90 -*-
!!
!! Filename: case.f90
!! Description: Interface between core and case code.
!! Author: Paul Bartholomew
!! Maintainer:
!! Created: Tue Feb 19 15:27:49 2019 (+0000)
!! Version:
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Commentary:
!!
!! This file defines an interface between the xcompact3d core and the
!! case-specific code, e.g. initialisation etc.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Change Log:
!!
!! [2019-02-19] Making module private by default.
!! [2019-02-19] Created case.f90
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! LICENSE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!! Code:

MODULE case

  USE param
  USE decomp_2d
  USE variables

  USE user_sim
  USE tgv
  USE cyl
  USE hill
  USE dbg_schemes
  USE channel
  USE mixlayer
  USE jet
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

    ELSEIF (itype.EQ.itype_hill) THEN

       CALL  init_hill (ux1,uy1,uz1,ep1,phi1)

    ELSEIF (itype.EQ.itype_cyl) THEN

       CALL init_cyl (ux1, uy1, uz1, phi1)

    ELSEIF (itype.EQ.itype_dbg) THEN

       CALL init_dbg (ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_mixlayer) THEN

       CALL init_mixlayer(rho1, ux1, uy1, uz1)

    ELSEIF (itype.EQ.itype_jet) THEN

       CALL init_jet(rho1, ux1, uy1, uz1, ep1, phi1)

    ELSEIF (itype.EQ.itype_tbl) THEN

       CALL init_tbl (ux1, uy1, uz1, ep1, phi1)

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

    ELSEIF (itype.EQ.itype_hill) THEN

       CALL boundary_conditions_hill (ux,uy,uz,phi,ep)

    ELSEIF (itype.EQ.itype_cyl) THEN

       CALL boundary_conditions_cyl (ux, uy, uz, phi)

    ELSEIF (itype.EQ.itype_dbg) THEN

       CALL boundary_conditions_dbg (ux, uy, uz, phi)

    ELSEIF (itype.EQ.itype_jet) THEN

       CALL boundary_conditions_jet (rho,ux,uy,uz,phi)

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

    ELSEIF (itype.EQ.itype_hill) THEN

       CALL postprocess_hill(ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_cyl) THEN

       CALL postprocess_cyl (ux, uy, uz, ep)

    ELSEIF (itype.EQ.itype_dbg) THEN

       CALL postprocess_dbg (ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_jet) THEN

       CALL postprocess_jet (ux, uy, uz, phi, ep)

    ELSEIF (itype.EQ.itype_tbl) THEN

      CALL postprocess_tbl (ux, uy, uz, pp, phi, ep)

    ENDIF

     if(iforces) then
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

    ELSEIF (itype.EQ.itype_jet) THEN

       CALL momentum_forcing_jet(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

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
