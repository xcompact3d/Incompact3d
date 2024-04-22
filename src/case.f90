!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module case

  use param
  use decomp_2d_constants
  use variables

  use case_generic
  use tgv
  use cyl
  use hill
  use channel
  use mixlayer
  use gravitycur
  use tbl
  use abl
  use uniform
  use sandbox
  use cavity
  use pipe

  use var, only : nzmsize

  implicit none

  logical :: case_visu_init = .false.
  
  private ! All functions/subroutines private by default
  public :: init, boundary_conditions, &
            momentum_forcing, scalar_forcing, set_fluid_properties, &
            test_flow, preprocessing, postprocessing, visu_case, visu_case_init

contains
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

    if (itype.eq.itype_generic) then

       call init_generic (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_gravitycur) then

       call init_gravitycur(rho1, ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_tgv) then

       call init_tgv (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_channel) then

       call init_channel (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_hill) then

       call  init_hill (ux1,uy1,uz1,ep1,phi1)

    elseif (itype.eq.itype_cyl) then

       call init_cyl (ux1, uy1, uz1, phi1)

    elseif (itype.eq.itype_mixlayer) then

       call init_mixlayer(rho1, ux1, uy1, uz1)

    elseif (itype.eq.itype_tbl) then

       call init_tbl (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_abl) then

       call init_abl (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_uniform) then

       call init_uniform (ux1, uy1, uz1, ep1, phi1)

    elseif (itype.EQ.itype_sandbox) THEN
   
       call init_sandbox (ux1, uy1, uz1, ep1, phi1, 0)

    elseif (itype.eq.itype_cavity) then

       call init_cavity(ux1, uy1, uz1, ep1, phi1)

    elseif (itype.eq.itype_pipe) then

       call init_pipe(ux1, uy1, uz1, ep1, phi1)

    else
  
         if (nrank.eq.0) then
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
    
    if (itype == itype_generic) then

       call boundary_conditions_generic (ux,uy,uz,phi,ep)

    elseif (itype.eq.itype_gravitycur) then

       call boundary_conditions_gravitycur(rho, phi)

    elseif (itype.eq.itype_tgv) then

       call boundary_conditions_tgv (ux, uy, uz, phi)

    elseif (itype.eq.itype_channel) then

       call boundary_conditions_channel (ux, uy, uz, phi)

    elseif (itype.eq.itype_hill) then

       call boundary_conditions_hill (ux,uy,uz,phi,ep)

    elseif (itype.eq.itype_cyl) then

       call boundary_conditions_cyl (ux, uy, uz, phi)

    elseif (itype.eq.itype_tbl) then

       call boundary_conditions_tbl (ux, uy, uz, phi)

    elseif (itype.eq.itype_abl) then

       call boundary_conditions_abl (ux, uy, uz, phi)

    elseif (itype.eq.itype_uniform) then

       call boundary_conditions_uniform (ux, uy, uz, phi)

    elseif (itype.EQ.itype_sandbox) THEN
   
       call boundary_conditions_sandbox (ux, uy, uz, phi)

    elseif (itype.eq.itype_cavity) then

       call boundary_conditions_cavity(ux, uy, uz, phi)

    elseif (itype.eq.itype_pipe) then

       call boundary_conditions_pipe (ux, uy, uz, phi)

    endif

  end subroutine boundary_conditions
  !##################################################################
  !##################################################################
  subroutine preprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d, only : xsize, ph1
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

    use decomp_2d, only : xsize, ph1
    use var, only : nzmsize, numscalar, nrhotime, npress, abl_T

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    integer :: j

    ! Recover temperature when decomposed (pressure to be recovered externally)
    if (itype.eq.itype_abl.and.ibuoyancy.eq.1) then
      do j=1,xsize(2) 
        abl_T(:,j,:,1) = phi1(:,j,:,1) + Tstat(j,1)
      enddo
      call run_postprocessing(rho1, ux1, uy1, uz1, pp3, abl_T, ep1)
    else
      call run_postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)
    endif

  end subroutine postprocessing
  !##################################################################
  !##################################################################
  subroutine run_postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d, only : xsize, ph1
    use visu, only  : write_snapshot, end_snapshot
    use stats, only : overall_statistic

    use var, only : nzmsize
    use var, only : itime
    use var, only : numscalar, nrhotime, npress

    use turbine, only : turbine_output
    use probes, only : write_probes

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    integer :: num

    if ((ivisu.ne.0).and.(mod(itime, ioutput).eq.0)) then
       call write_snapshot(rho1, ux1, uy1, uz1, pp3, phi1, ep1, itime, num)

       ! XXX: Ultimate goal for ADIOS2 is to pass do all postproc online - do we need this?
       !      Currently, needs some way to "register" variables for IO
       call visu_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1, num)

       call end_snapshot(itime, num)
    end if

    call postprocess_case(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    call overall_statistic(ux1, uy1, uz1, phi1, pp3, ep1)

    if (iturbine.ne.0) then 
      call turbine_output()
    endif

    call write_probes(ux1, uy1, uz1, pp3, phi1)

  end subroutine run_postprocessing
  !##################################################################
  !##################################################################
  subroutine postprocess_case(rho,ux,uy,uz,pp,phi,ep)

    use forces
    use var, only : nzmsize
    use param, only : npress

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp

    if (itype.eq.itype_generic) then

       call postprocess_generic (ux, uy, uz, phi, ep)

    elseif (itype.eq.itype_gravitycur) then

       call postprocess_gravitycur(rho, ux, uy, uz, phi, ep)

    elseif (itype.eq.itype_tgv) then

       call postprocess_tgv (ux, uy, uz, phi, ep)

    elseif (itype.eq.itype_channel) then

       call postprocess_channel (ux, uy, uz, pp, phi, ep)

    elseif (itype.eq.itype_hill) then

       call postprocess_hill(ux, uy, uz, pp, phi, ep)

    elseif (itype.eq.itype_cyl) then

       call postprocess_cyl (ux, uy, uz, ep)

    elseif (itype.eq.itype_tbl) then

       call postprocess_tbl (ux, uy, uz, ep)

    elseif (itype.eq.itype_abl) then

       call postprocess_abl (ux, uy, uz, ep)

    elseif (itype.eq.itype_uniform) then

       call postprocess_uniform (ux, uy, uz, ep)

    elseif (itype.EQ.itype_sandbox) THEN
   
       call postprocess_sandbox (ux, uy, uz, phi, ep)

    elseif (itype.eq.itype_cavity) then

       call postprocess_cavity(ux, uy, uz, phi)

    elseif (itype.eq.itype_pipe) then

       call postprocess_pipe(ux, uy, uz, pp, phi, ep)

    endif

    if (iforces.eq.1) then
       call force(ux,uy,ep)
       call restart_forces(1)
    endif

  end subroutine postprocess_case
  !##################################################################
  !!
  !!  SUBROUTINE: visu_case_init
  !!      AUTHOR: PB
  !! DESCRIPTION: Initialise case-specific visualization
  !!
  !##################################################################
  subroutine visu_case_init

    implicit none
    
    if (itype .eq. itype_tgv) then

       call visu_tgv_init(case_visu_init)

    else if (itype .eq. itype_channel) then

       call visu_channel_init(case_visu_init)

    else if (itype .eq. itype_hill) then

       call visu_hill_init(case_visu_init)

    else if (itype .eq. itype_cyl) then

       call visu_cyl_init(case_visu_init)

    else if (itype .eq. itype_tbl) then

       call visu_tbl_init(case_visu_init)

    else if (itype .eq. itype_gravitycur) then

       call visu_gravitycur_init(case_visu_init)

    else if (itype .eq. itype_uniform) then

       call visu_uniform_init(case_visu_init)      

    end if
    
  end subroutine visu_case_init
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
    integer, intent(in) :: num

    logical :: called_visu = .false.
    
    if (itype.eq.itype_generic) then

       call visu_generic(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.
       
    elseif (itype.eq.itype_tgv) then

       call visu_tgv(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.

    elseif (itype.eq.itype_channel) then

       call visu_channel(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.

    elseif (itype.eq.itype_hill) then

       call visu_hill(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.

    elseif (itype.eq.itype_cyl) then

       call visu_cyl(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.

    elseif (itype.eq.itype_tbl) then

       call visu_tbl(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.
       
    elseif (itype.eq.itype_uniform) then

       call visu_uniform(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.

    endif

    if (called_visu .and. (.not. case_visu_init)) then

       print *, "ERROR: tried to run case-specific visu without initialisation!"
       print *, "       See the TGV case initialisation for example."
       STOP
       
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

    use mhd, only: mhd_active,momentum_forcing_mhd

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    if (itype.eq.itype_channel) then

       call momentum_forcing_channel(dux1, duy1, duz1, ux1, uy1, uz1)

    elseif (itype.eq.itype_abl) then

       call momentum_forcing_abl(dux1, duy1, duz1, ux1, uy1, uz1, phi1)

    endif

    if(mhd_active) then
      call momentum_forcing_mhd(dux1(:,:,:,1),duy1(:,:,:,1),duz1(:,:,:,1),ux1,uy1,uz1)
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
  subroutine scalar_forcing(dphi1, rho1, ux1, uy1, uz1, phi1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, phi1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphi1

    if (itype.eq.itype_abl) then

       call scalar_forcing_abl(uy1, dphi1, phi1)

    endif

  end subroutine scalar_forcing
  !##################################################################
  !##################################################################
  subroutine set_fluid_properties(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

    if (itype.eq.itype_gravitycur) then

       if (ilmn) then 
          call set_fluid_properties_gravitycur(rho1, mu1)
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
  !##################################################################
end module case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
