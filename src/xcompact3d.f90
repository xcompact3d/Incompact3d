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

program xcompact3d

  use var
  use case

  use param, only : catching_signal
  use transeq, only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
       calc_divu_constraint, solve_poisson, cor_vel
  use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
  use turbine, only : compute_turbines
  use ibm_param
  use ibm, only : body
  use genepsi, only : genepsi3d

  implicit none

  call init_xcompact3d()

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     call simu_stats(2)

     if (iturbine /= 0) call compute_turbines(ux1, uy1, uz1)

     if (iin == 3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype == itype_abl.or.iturbine /= 0).and.(ifilter /= 0).and.(ilesmod /= 0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

        if (imove == 1) then ! update epsi for moving objects
          if ((iibm == 2).or.(iibm == 3)) then
            call genepsi3d(ep1)
          else if (iibm == 1) then
            call body(ux1,uy1,uz1,ep1)
          endif
        endif

        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
           call velocity_to_momentum(rho1,ux1,uy1,uz1)
        endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

     enddo !! End sub timesteps

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

     if (catching_signal /= 0) call catch_signal()

  enddo !! End time loop

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use forces

  use var

  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
                    restart, simu_stats, compute_cfldiff

  use param, only : ilesmod, jles,itype
  use param, only : irestart, catching_signal

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : nstat, nvisu, nprobe

  use les, only: init_explicit_les
  use turbine, only: init_turbines

  use visu, only : visu_init

  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body

  use probes, only : init_probes

  implicit none

  external catch_sigusr1, catch_sigusr2

  ! This should be valid for ARM and X86. See "man 7 signal"
  integer, parameter :: SIGUSR1 = 10, SIGUSR2 = 12
  integer :: code

  integer :: nargin, FNLength, status
  character(len=80) :: InputFN

  !! Initialise MPI
  call MPI_INIT(code)
  if (code /= 0) call decomp_2d_abort(code, "MPI_INIT")
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,code)
  if (code /= 0) call decomp_2d_abort(code, "MPI_COMM_RANK")
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)
  if (code /= 0) call decomp_2d_abort(code, "MPI_COMM_SIZE")

  ! Install signal handler
  call signal(SIGUSR1, catch_sigusr1)
  call signal(SIGUSR2, catch_sigusr2)
  catching_signal = 0

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) print*, 'Xcompact3d is run with the default file -->', InputFN
  elseif (nargin >= 1) then
     if (nrank==0) print*, 'Xcompact3d is run with the provided file -->', InputFN

     call get_command_argument(1,InputFN,FNLength,status)
     if (status /= 0) then
        if (nrank == 0) print*, 'InputFN is too small for the given input file'
        call decomp_2d_abort(status, "get_command_argument")
     endif
     if (nrank==0) print*, 'Xcompact3d is run with the provided file -->', InputFN
  endif

#ifdef ADIOS2
  if (nrank  ==  0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod /= 0) then
     if (jles > 0)  call init_explicit_les()
  endif

  if ((iibm == 2).or.(iibm == 3)) then
     call genepsi3d(ep1)
  else if (iibm == 1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (iforces == 1) then
     call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif

  !####################################################################
  ! initialise visu
  if (ivisu /= 0) call visu_init()
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,0)
  endif

  if ((iibm == 2).or.(iibm == 3)) then
     call genepsi3d(ep1)
  else if (iibm == 1) then
     call body(ux1,uy1,uz1,ep1)
  endif

  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()

  if (iturbine /= 0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank == 0)then
        open(42,file='time_evol.dat',form='formatted')
     endif
  endif
  if (itype==5) then
     if(nrank == 0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_poisson, only : decomp_2d_poisson_finalize
  use tools, only : simu_stats
  use var, only : finalize_variables
  use forces, only : iforces, finalize_forces
  use param, only : iimplicit, itype
  use ydiff_implicit, only : finalize_implicit
  use case, only : finalize_case
  use probes, only : finalize_probes
  use visu, only : finalize_visu

  implicit none

  integer :: code
  
  if (itype==2) then
     if(nrank == 0)then
        close(42)
     endif
  endif
  if (itype==5) then
     if(nrank == 0)then
        close(38)
     endif
  endif
  
  call simu_stats(4)
  if (iforces == 1) call finalize_forces()
  if (iimplicit /= 0) call finalize_implicit()
  call finalize_case()
  call finalize_probes()
  call finalize_visu()
  call finalize_variables()
  call decomp_info_finalize(ph1)
  call decomp_info_finalize(ph4)
  call decomp_info_finalize(ph2)
  call decomp_info_finalize(ph3)
  call decomp_info_finalize(phG)
  call decomp_2d_finalize()
  call decomp_2d_poisson_finalize()
  CALL MPI_FINALIZE(code)
  if (code /= 0) call decomp_2d_abort(code, "MPI_FINALIZE")

endsubroutine finalise_xcompact3d
!########################################################################
!########################################################################
subroutine catch_sigusr1

  use decomp_2d, only : nrank
  use param, only : catching_signal

  implicit none

  if (nrank == 0) print *, "Catching SIGUSR1"
  catching_signal = 1

end subroutine catch_sigusr1
!
subroutine catch_sigusr2

  use decomp_2d, only : nrank
  use param, only : catching_signal

  implicit none

  if (nrank == 0) print *, "Catching SIGUSR2"
  catching_signal = 2

end subroutine catch_sigusr2
!
subroutine catch_signal

  use MPI
  use decomp_2d, only : decomp_2d_abort
  use param, only : catching_signal, itime, ilast, icheckpoint

  implicit none

  integer :: code

  if (catching_signal == 1) then

    ! Stop at the end of the next iteration
    ilast = itime + 1
    call MPI_ALLREDUCE(MPI_IN_PLACE, ilast, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

    return

  elseif (catching_signal == 2) then

    ! Stop at the end of the next iteration
    ilast = itime + 1
    call MPI_ALLREDUCE(MPI_IN_PLACE, ilast, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, code)
    if (code /= 0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

    ! Write checkpoint at the end of the next iteration
    icheckpoint = ilast

    return

  endif

end subroutine catch_signal
