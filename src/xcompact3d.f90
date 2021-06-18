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

     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype.eq.itype_abl.or.iturbine.ne.0).and.(ifilter.ne.0).and.(ilesmod.ne.0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

if (imove.eq.1) then ! update epsi for moving objects
  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call body(ux1,uy1,uz1,ep1)
  endif
endif
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

        !! XXX N.B. from this point, X-pencil velocity arrays contain momentum.
        call velocity_to_momentum(rho1,ux1,uy1,uz1)

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        call momentum_to_velocity(rho1,ux1,uy1,uz1)
        !! XXX N.B. from this point, X-pencil velocity arrays contain velocity.

        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

     enddo !! End sub timesteps

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)

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
       restart, &
       simu_stats, compute_cfldiff

  use param, only : ilesmod, jles,itype
  use param, only : irestart

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

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase

  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) print*, 'Xcompact3d is run with the default file -->', InputFN
  elseif (nargin.ge.1) then
     if (nrank==0) print*, 'Program is run with the provided file -->', InputFN

     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
  endif

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

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (iforces.eq.1) then
     call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif

  !####################################################################
  ! initialise visu
  if (ivisu.ne.0) call visu_init()
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,0)
!     ux1(:,:,:)=ux1(:,:,:)-0.5
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if ((iibm.eq.1).or.(iibm.eq.3)) then
     call body(ux1,uy1,uz1,ep1)
  endif

  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes(ep1)

  if (iturbine.ne.0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank.eq.0)then
        open(42,file='time_evol.dat',form='formatted')
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d

  use tools, only : simu_stats
  use param, only : itype

  implicit none

  integer :: ierr
  
  if (itype==2) then
     if(nrank.eq.0)then
        close(42)
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        close(38)
     endif
  endif
  
  call simu_stats(4)
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d
