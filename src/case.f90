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
!! This file defines an interface between the incompact3d core and the
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
  
  USE tgv
  USE cyl
  USE hill
  USE dbg_schemes
  USE channel
  USE mixlayer

  IMPLICIT NONE

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init, boundary_conditions, postprocessing

CONTAINS

  SUBROUTINE init (rho1, ux1, uy1, uz1, ep1, phi1, drho1, dux1, duy1, duz1, phis1, phiss1)

    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar) :: phi1,phis1,phiss1
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1, drho1

    !! Default density and pressure0 to one
    pressure0 = one
    rho1(:,:,:,:) = one

    IF (itype.EQ.itype_lockexch) THEN

       IF (nrank.EQ.0) THEN
          PRINT *, "Lock-exchange case not modernised yet!"
          STOP
       ENDIF
       
    ELSEIF (itype.EQ.itype_tgv) THEN
       
       CALL init_tgv (ux1, uy1, uz1, ep1, phi1, dux1, duy1, duz1, phis1, phiss1)
       
    ELSEIF (itype.EQ.itype_channel) THEN
       
       CALL init_channel (ux1, uy1, uz1, ep1, phi1, dux1, duy1, duz1, phis1, phiss1)
       
    ELSEIF (itype.EQ.itype_hill) THEN

       CALL  init_hill (ux1,uy1,uz1,ep1,phi1,dux1,duy1,duz1,phis1,phiss1)
          
    ELSEIF (itype.EQ.itype_cyl) THEN
       
       CALL init_cyl (ux1, uy1, uz1, ep1, phi1, dux1, duy1, duz1, phis1, phiss1)
       
    ELSEIF (itype.EQ.itype_dbg) THEN
       
       CALL init_dbg (ux1, uy1, uz1, ep1, phi1, dux1, duy1, duz1, phis1, phiss1)

    ELSEIF (itype.EQ.itype_mixlayer) THEN

       CALL init_mixlayer(rho1, ux1, uy1, uz1, drho1, dux1, duy1, duz1)
       
    ENDIF

  END SUBROUTINE init

  SUBROUTINE boundary_conditions (ux,uy,uz,phi,ep)
    
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    IF (itype.EQ.itype_lockexch) THEN

       IF (nrank.EQ.0) THEN
          PRINT *, "Lock-exchange case not modernised yet!"
          STOP
       ENDIF
       
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
       
    ENDIF
    
  END SUBROUTINE boundary_conditions

  SUBROUTINE postprocessing(ux,uy,uz,phi,ep)
    
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    REAL(mytype),DIMENSION(xsize(1),xsize(2),xsize(3)) :: ep

    IF (itype.EQ.itype_lockexch) THEN

       IF (nrank.EQ.0) THEN
          PRINT *, "Lock-exchange case not modernised yet!"
          STOP
       ENDIF
       
    ELSEIF (itype.EQ.itype_tgv) THEN
       
       CALL postprocessing_tgv (ux, uy, uz, phi, ep)
       
    ELSEIF (itype.EQ.itype_channel) THEN
       
       CALL postprocessing_channel (ux, uy, uz, phi, ep)
       
    ELSEIF (itype.EQ.itype_hill) THEN
       
       CALL postprocessing_hill(ux, uy, uz, phi, ep)
       
    ELSEIF (itype.EQ.itype_cyl) THEN
       
       CALL postprocessing_cyl (ux, uy, uz, phi, ep)
       
    ELSEIF (itype.EQ.itype_dbg) THEN
       
       CALL postprocessing_dbg (ux, uy, uz, phi, ep)
       
    ENDIF
  END SUBROUTINE postprocessing
  
END MODULE case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
