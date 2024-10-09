!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module fast_projection_method

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d

  implicit none

  private
  public :: init_fast_projection_method, solve_pressure_fast_projection

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: init_fast_projection_method
  !! DESCRIPTION: Initialize module specific parameters and variables
  !!      INPUTS: 
  !!     OUTPUTS: 
  !!      AUTHOR: Kay Schäfer, Carlo de Michele
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fast_projection_method()

    use param
    use variables

    implicit none

    integer :: i, j, k
    
    !! TODO: fill routine

  end subroutine init_fast_projection_method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: solve_pressure_fast_projection
  !! DESCRIPTION: Manages the computation of the fast-projection method
  !!              depending on the Runge-Kutta sub-time step.
  !!      INPUTS: 
  !!     OUTPUTS: 
  !!      AUTHOR: Kay Schäfer, Carlo de Michele
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine solve_pressure_fast_projection(pp3, px1, py1, pz1, rho1, ux1, uy1, uz1, ep1, drho1, divu3)

    use decomp_2d_poisson, ONLY : poisson
    use var, ONLY : nzmsize
    use param, ONLY : itr, ntime, nrhotime, npress
    use param, ONLY : zero, one, two
  
    use navier, only : solve_poisson, gradp, calc_divu_constraint

    implicit none

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime), INTENT(IN) :: drho1
    
    REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(INOUT) :: divu3

    !! Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), INTENT(INOUT) :: pp3
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    ! if (nrank==0) then
    !   write(*,*) '==========================================================='
    !   write(*,"(' itr           : ',I9)") itr
    !   write(*,*) '==========================================================='
    ! end if

    !! Check for sub-time step in Runge-Kutta method 3rd order 
    if ((itr.eq.one).or.(itr.eq.two)) then 
      !! Extrapolate pressure field from previous two pressure instances
      call extrapolate_time_pressure(pp3)
      !! Alternatively one can manage the if-conditions for itr in extrapolate_time_pressure here
      !! and use the routine extrapolate_time_pressure(pp3,dt1,dt2) only to extrapolate for the 
      !! provided pressure fields and time increments 


      !! Compute pressure gradients for pressure correction step 
      CALL gradp(px1,py1,pz1,pp3(:,:,:,1))

    else 
      !! Last RK3 step: solve poisson to get divergence-free velocity field 
      call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
    end if 

  end subroutine solve_pressure_fast_projection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: extrapolate_time_pressure
  !! DESCRIPTION: Extrapolates the pressure by previous pressure fields in time
  !!      INPUTS: 
  !!     OUTPUTS: 
  !!      AUTHOR: Kay Schäfer, Carlo de Michele
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine extrapolate_time_pressure(pp3)

    use param, ONLY : zero, one, two, npress, dt, itr
    use variables
    use var, ONLY: ta1, nzmsize

    implicit none

    integer :: it, i, j, k
    real(mytype) :: dt_substep1, dt_substep2
    
    !! In/Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), INTENT(INOUT) :: pp3

    !! TODO: 
    !!    - determine dt_substep1/2 for each time step, since we want to have
    !!      the possibility for a flexible time adjustment
    !!    - maybe roll out for-loops for computation
    !!

    !! if itr==1 then extrapolate from p^(n-1) and p^n to first RK3 substep 
    !! Compute 1/dt before and save it to scalar value!
    if (itr.eq.one) then 
      ta1 = ( pp3(:,:,:,1) - pp3(:,:,:,2) ) / (dt) * dt_substep1 + pp3(:,:,:,1)

    elseif (itr.eq.two) then
      !! if itr==2 then extrapolate from p^n and first RK3 substep  p^n_(k=1)  to second RK3 substep
      ta1 = ( pp3(:,:,:,1) - pp3(:,:,:,2) ) / (dt_substep1) * dt_substep2 + pp3(:,:,:,1)
    end if 

    !! Flip pressure data two avoid a third array for previous pressure field
    pp3(:,:,:,2) = pp3(:,:,:,1)
    pp3(:,:,:,1) = ta1(:,:,:)


  end subroutine extrapolate_time_pressure

end module fast_projection_method
