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

!########################################################################
subroutine boot_xcompact3d()

  use MPI
  use decomp_2d, only : nrank, nproc, decomp_2d_abort
  use param, only : catching_signal

  implicit none

  ! This should be valid for ARM and X86. See "man 7 signal"
  external catch_sigusr1, catch_sigusr2
  integer, parameter :: SIGUSR1 = 10, SIGUSR2 = 12
  integer :: code

  integer :: nargin, FNLength, status
  character(len=80) :: InputFN

  !! Initialise MPI
  call MPI_INIT(code)
  if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_INIT")
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,code)
  if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_RANK")
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)
  if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_SIZE")

  ! Install signal handler
  call signal(SIGUSR1, catch_sigusr1)
  call signal(SIGUSR2, catch_sigusr2)
  catching_signal = 0

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'Xcompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     if (status /= 0) then
        if (nrank == 0) print*, 'InputFN is too small for the given input file'
        call decomp_2d_abort(__FILE__, __LINE__, status, "get_command_argument")
     endif
  endif

#ifdef ADIOS2
  if (nrank == 0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call init_xcompact3d(inputFN)

endsubroutine boot_xcompact3d
!########################################################################

!########################################################################
subroutine init_xcompact3d(InputFN)

  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_init
  use decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use sandbox, only : init_sandbox
  use forces
  use var
  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
                    restart, simu_stats, compute_cfldiff, &
                    init_inflow_outflow
  use param, only : ilesmod, jles,itype, &
                    irestart, catching_signal
  use variables, only : nx, ny, nz, nxm, nym, nzm, &
                        p_row, p_col, &
                        nstat, nvisu, nprobe
  use les, only: init_explicit_les
  use turbine, only: init_turbines
  use visu, only : visu_init, visu_ready
  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body
  use probes, only : init_probes

  implicit none

  character(len=*), intent(in) :: InputFN
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_io_init()
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
  if (ivisu /= 0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  call boot(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     endif
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif

  if ((ioutflow == 1).or.(iin == 3)) then
     call init_inflow_outflow()
  end if

  if ((iibm == 2).or.(iibm == 3)) then
     call genepsi3d(ep1)
  else if (iibm == 1) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()

  if (iturbine /= 0) call init_turbines(ux1, uy1, uz1)

endsubroutine init_xcompact3d
!########################################################################

!########################################################################
subroutine finalise_xcompact3d(flag)

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise
  use decomp_2d_poisson, only : decomp_2d_poisson_finalize
  use tools, only : simu_stats
  use var, only : finalize_variables
  use forces, only : iforces, finalize_forces
  use param, only : iimplicit, itype
  use ydiff_implicit, only : finalize_implicit
  use case, only : finalize_case
  use probes, only : finalize_probes
  use visu, only : visu_finalise

  implicit none

  ! Set to false to skip MPI_FINALIZE (unit tests)
  logical, intent(in) :: flag

  integer :: code

  call simu_stats(4)
  if (iforces == 1) call finalize_forces()
  if (iimplicit /= 0) call finalize_implicit()
  call finalize_case()
  call finalize_probes()
  call visu_finalise()
  call finalize_variables()
  call decomp_info_finalize(ph1)
  call decomp_info_finalize(ph4)
  call decomp_info_finalize(ph2)
  call decomp_info_finalize(ph3)
  call decomp_info_finalize(phG)
  call decomp_2d_io_finalise()
  call decomp_2d_finalize()
  call decomp_2d_poisson_finalize()
  if (flag) then
    CALL MPI_FINALIZE(code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_FINALIZE")
  endif

endsubroutine finalise_xcompact3d
!########################################################################

subroutine check_transients()

  use decomp_2d, only : mytype

  use var
  use tools, only : avg3d

  implicit none

  real(mytype) avg_param

  avg_param = zero
  call avg3d (dux1, avg_param)
  if (nrank == 0) write(*,*)'## Main dux1 ', avg_param
  avg_param = zero
  call avg3d (duy1, avg_param)
  if (nrank == 0) write(*,*)'## Main duy1 ', avg_param
  avg_param = zero
  call avg3d (duz1, avg_param)
  if (nrank == 0) write(*,*)'## Main duz1 ', avg_param

end subroutine check_transients
!
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
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE")

    catching_signal = 0
    return

  elseif (catching_signal == 2) then

    ! Stop at the end of the next iteration
    ilast = itime + 1
    call MPI_ALLREDUCE(MPI_IN_PLACE, ilast, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE")

    ! Write checkpoint at the end of the next iteration
    icheckpoint = ilast

    catching_signal = 0
    return

  endif

end subroutine catch_signal
