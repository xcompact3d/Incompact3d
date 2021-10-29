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

module case

  use param
  use decomp_2d
  use variables

  use user_sim
  use tgv
  use cyl
  use hill
  use dbg_schemes
  use channel
  use mixlayer
  use jet
  use lockexch
  use tbl
  use abl
  use uniform
  use sandbox
  use cavity

  use var, only : nzmsize

  implicit none

  private ! All functions/subroutines private by default
  public :: boot, init, boundary_conditions, &
            momentum_forcing, scalar_forcing, set_fluid_properties, &
            test_flow, preprocessing, postprocessing, finalize_case, &
            visu_case

contains
  !##################################################################
  subroutine boot(rho1, ux1, uy1, uz1, ep1, phi1, drho1, dux1, duy1, duz1, dphi1, &
       pp3, px1, py1, pz1)

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1,drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1

    if (itype == itype_tgv) then

       call boot_tgv()

    elseif (itype == itype_channel) then

       call boot_channel()

    elseif (itype == itype_cyl) then

       call boot_cyl()

    endif

  end subroutine boot
  !##################################################################
  !##################################################################
  subroutine init (rho1, ux1, uy1, uz1, ep1, phi1, drho1, dux1, duy1, duz1, dphi1, &
       pp3, px1, py1, pz1)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1,drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1

    INTEGER :: it, is

    !! Zero out the pressure field
    pp3(:,:,:,1) = zero
    px1(:,:,:) = zero
    py1(:,:,:) = zero
    pz1(:,:,:) = zero

    !! Default density and pressure0 to one
    pressure0 = one
    rho1(:,:,:,:) = one

    if (itype == itype_user) then

       call init_user (ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_lockexch) then

       call init_lockexch(rho1, ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_tgv) then

       call init_tgv (ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_channel) then

       call init_channel (ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_hill) then

       call  init_hill (ux1,uy1,uz1,ep1,phi1)

    elseif (itype == itype_cyl) then

       call init_cyl (ux1, uy1, uz1, phi1)

    elseif (itype == itype_dbg) then

       call init_dbg (ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_mixlayer) then

       call init_mixlayer(rho1, ux1, uy1, uz1)

    elseif (itype == itype_jet) then

       call init_jet(rho1, ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_tbl) then

       call init_tbl (ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_abl) then

       call init_abl (ux1, uy1, uz1, ep1, phi1)

    elseif (itype == itype_uniform) then

       call init_uniform (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.EQ.itype_sandbox) THEN
   
       call init_sandbox (ux1, uy1, uz1, ep1, phi1, 0)

    elseif (itype == itype_cavity) then

       call init_cavity (ux1, uy1, uz1, ep1, phi1)

    else
  
         if (nrank == 0) then
            print *, "ERROR: Unknown itype: ", itype
            STOP
         endif

    endif

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

  end subroutine init
  !##################################################################
  !##################################################################
  subroutine boundary_conditions (rho,ux,uy,uz,phi,ep)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho

    if (itype == itype_user) then

       call boundary_conditions_user (ux,uy,uz,phi,ep)

    elseif (itype == itype_lockexch) then

       call boundary_conditions_lockexch(rho, phi)

    elseif (itype == itype_tgv) then

       call boundary_conditions_tgv (ux, uy, uz, phi)

    elseif (itype == itype_channel) then

       call boundary_conditions_channel (ux, uy, uz, phi)

    elseif (itype == itype_hill) then

       call boundary_conditions_hill (ux,uy,uz,phi,ep)

    elseif (itype == itype_cyl) then

       call boundary_conditions_cyl (ux, uy, uz, phi)

    elseif (itype == itype_dbg) then

       call boundary_conditions_dbg (ux, uy, uz, phi)

    elseif (itype == itype_jet) then

       call boundary_conditions_jet (rho,ux,uy,uz,phi)

    elseif (itype == itype_tbl) then

       call boundary_conditions_tbl (ux, uy, uz, phi)

    elseif (itype == itype_abl) then

       call boundary_conditions_abl (ux, uy, uz, phi)

    elseif (itype == itype_uniform) then

       call boundary_conditions_uniform (ux, uy, uz, phi)

    elseif (itype.EQ.itype_sandbox) THEN
   
       call boundary_conditions_sandbox (ux, uy, uz, phi)

    elseif (itype == itype_cavity) then

       call boundary_conditions_cavity (ux, uy, uz, phi)

    endif

  end subroutine boundary_conditions
  !##################################################################
  !##################################################################
  subroutine preprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d, only : mytype, xsize, ph1
    use visu, only  : write_snapshot
    use stats, only : overall_statistic

    use var, only : nzmsize
    use var, only : itime
    use var, only : numscalar, nrhotime, npress

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    !call write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime)
    !call postprocess_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1)
    !call overall_statistic(ux1, uy1, uz1, phi1, pp3, ep1)

  end subroutine preprocessing
  !##################################################################
  !##################################################################
  subroutine postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d, only : mytype, xsize, ph1
    use visu, only  : write_snapshot, end_snapshot
    use stats, only : overall_statistic

    use var, only : nzmsize
    use var, only : itime
    use var, only : numscalar, nrhotime, npress
    use var, only : T_tmp

    use turbine, only : turbine_output
    use probes, only : write_probes

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(inout) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    integer :: i, j, k
    character(len=32) :: num

    ! Recover temperature when decomposed (pressure to be recovered externally)
    if (itype == itype_abl.and.ibuoyancy == 1) then
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          do i = 1, xsize(1)
            T_tmp(i,j,k,1) = phi1(i,j,k,1)
            phi1(i,j,k,1) = phi1(i,j,k,1) + Tstat(j,1)
          enddo
        enddo
      enddo
    endif

    if ((ivisu /= 0).and.(mod(itime, ioutput) == 0)) then
      call write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime, num)
#ifndef ADIOS2
       ! XXX: Ultimate goal for ADIOS2 is to pass do all postproc online - do we need this?
       !      Currently, needs some way to "register" variables for IO
      call visu_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1, num)
#endif
      call end_snapshot(itime, num)
    end if

    call postprocess_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    call overall_statistic(ux1, uy1, uz1, phi1, pp3, ep1)

    if (iturbine /= 0) then 
      call turbine_output()
    endif

    call write_probes(ux1, uy1, uz1, pp3, phi1)

    if (itype == itype_abl.and.ibuoyancy == 1) then
      phi1(:,:,:,1) = T_tmp(:,:,:,1)
    endif

  end subroutine postprocessing
  !##################################################################
  !##################################################################
  subroutine postprocess_case(rho,ux,uy,uz,pp,phi,ep)

    use forces, only : iforces, force, restart_forces
    use var, only : nzmsize
    use param, only : npress

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp

    if (itype == itype_user) then

       call postprocess_user (ux, uy, uz, phi, ep)

    elseif (itype == itype_lockexch) then

       call postprocess_lockexch(rho, ux, uy, uz, phi, ep)

    elseif (itype == itype_tgv) then

       call postprocess_tgv (ux, uy, uz, phi, ep)

    elseif (itype == itype_channel) then

       call postprocess_channel (ux, uy, uz, pp, phi, ep)

    elseif (itype == itype_hill) then

       call postprocess_hill(ux, uy, uz, phi, ep)

    elseif (itype == itype_cyl) then

       call postprocess_cyl (ux, uy, uz, ep)

    elseif (itype == itype_dbg) then

       call postprocess_dbg (ux, uy, uz, phi, ep)

    elseif (itype == itype_jet) then

       call postprocess_jet (ux, uy, uz, phi, ep)

    elseif (itype == itype_tbl) then

       call postprocess_tbl (ux, uy, uz, ep)

    elseif (itype == itype_abl) then

       call postprocess_abl (ux, uy, uz, ep)

    elseif (itype == itype_uniform) then

       call postprocess_uniform (ux, uy, uz, ep)

    elseif (itype.EQ.itype_sandbox) THEN
   
       call postprocess_sandbox (ux, uy, uz, phi, ep)

    elseif (itype == itype_cavity) then

       call postprocess_cavity(ux, uy, uz, pp, phi)

    endif

    if (iforces == 1) then
       if (itype == itype_cyl) call force(ux, uy, ep, cyl_iounit)
       call restart_forces(1)
    endif

  end subroutine postprocess_case
  !##################################################################
  !!
  !!  SUBROUTINE: visu_case
  !!      AUTHOR: CF
  !! DESCRIPTION: Call case-specific visualization
  !!
  !##################################################################
  subroutine visu_case(rho1,ux1,uy1,uz1,pp3,phi1,ep1,num)

    use var, only : nzmsize
    use param, only : npress

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

    if (itype == itype_user) then

       call visu_user(ux1, uy1, uz1, pp3, phi1, ep1, num)

    elseif (itype == itype_tgv) then

       call visu_tgv(ux1, uy1, uz1, pp3, phi1, ep1, num)

    elseif (itype == itype_channel) then

       call visu_channel(ux1, uy1, uz1, pp3, phi1, ep1, num)

    elseif (itype == itype_cyl) then

       call visu_cyl(ux1, uy1, uz1, pp3, phi1, ep1, num)

    elseif (itype == itype_tbl) then

       call visu_tbl(ux1, uy1, uz1, pp3, phi1, ep1, num)

    endif

  end subroutine visu_case
  !##################################################################
  !##################################################################
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Calls case-specific forcing functions for the
  !!              momentum equations.
  !!
  !##################################################################
  subroutine momentum_forcing(dux1, duy1, duz1, rho1, ux1, uy1, uz1, phi1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    if (itype == itype_channel) then

       call momentum_forcing_channel(dux1, duy1, duz1, ux1, uy1, uz1)

    elseif (itype == itype_jet) then

       call momentum_forcing_jet(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

    elseif (itype == itype_abl) then

       call momentum_forcing_abl(dux1, duy1, duz1, ux1, uy1, uz1, phi1)

    elseif (itype == itype_cavity) then

       call momentum_forcing_cavity(duy1, phi1)

    endif

  end subroutine momentum_forcing
  !##################################################################
  !##################################################################
  !!
  !!  SUBROUTINE: scalar_forcing
  !!      AUTHOR: Kay Schäfer
  !! DESCRIPTION: Calls case-specific forcing functions for the
  !!              scalar transport equations.
  !!
  !##################################################################
  subroutine scalar_forcing(dphi1, rho1, ux1, uy1, uz1, phi1, is)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, phi1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphi1
    integer, optional, intent(in) :: is

    if (itype == itype_channel) then

       call scalar_forcing_channel(dphi1, phi1, is)

    elseif (itype == itype_abl) then

       call scalar_forcing_abl(uy1, dphi1, phi1)

    endif

  end subroutine scalar_forcing
  !##################################################################
  !##################################################################
  subroutine set_fluid_properties(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    if (itype == itype_lockexch) then

       if (ilmn) then 
          call set_fluid_properties_lockexch(rho1, mu1)
       end if
       
    endif

  endsubroutine set_fluid_properties
  !##################################################################
  !##################################################################
  subroutine test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

    use decomp_2d
    use param

    use navier, only : divergence

    use var, only : numscalar, dv3
    use tools, only : test_speed_min_max, compute_cfl, test_scalar_min_max
    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1, ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime), intent(in) :: drho1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    if ((mod(itime,ilist)==0 .or. itime == ifirst .or. itime == ilast)) then
       call divergence(dv3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,2)
       call test_speed_min_max(ux1,uy1,uz1)
       call compute_cfl(ux1,uy1,uz1)
       if (iscalar==1) call test_scalar_min_max(phi1)
    endif

  end subroutine test_flow
  !##################################################################
  !!
  !!  SUBROUTINE: finalize_case
  !!      AUTHOR: Cédric Flageul
  !! DESCRIPTION: Calls case-specific functions to close units
  !!              and free memory
  !!
  !##################################################################
  subroutine finalize_case()

    implicit none

    if (itype == itype_tgv) then

       call finalize_tgv()

    elseif (itype == itype_channel) then

       call finalize_channel()

    elseif (itype == itype_cyl) then

       call finalize_cyl()

    elseif (itype == itype_cavity) then

       call finalize_cavity()

    endif

  end subroutine finalize_case
end module case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
